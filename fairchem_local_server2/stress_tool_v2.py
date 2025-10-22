"""
UMA WS stress utility — MIXED + ERROR DETAILS (plus)
- Everything from v3_mixed_errors, plus:
  * Synthetic error event when a session ends with NO first reply (phase='no_reply').
  * Records non-1000 clean closes as errors (phase='close').
  * Adds --open-timeout (seconds) to tune WebSocket handshake.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import random
import time
import traceback as tb
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple

import websockets
import websockets.exceptions as wse

from fairchem_local_server2 import session_pb2 as pb

# ---- tiny geometry seed -----------------------------------------------------


def _water():
    Z = [8, 1, 1]
    xyz = [
        [0.0000, 0.0000, 0.0000],
        [0.9572, 0.0000, 0.0000],
        [-0.239987, 0.926627, 0.0000],
    ]
    return Z, xyz


# ---- metrics ----------------------------------------------------------------


@dataclass
class Metrics:
    sessions_ok: int = 0
    sessions_err: int = 0
    messages_rcvd: int = 0
    first_reply_lat_ms_accum: float = 0.0
    first_reply_lat_count: int = 0
    frames_proto_ok: int = 0
    frames_proto_parse_err: int = 0
    frames_text: int = 0
    frames_json_in_bytes: int = 0
    frames_other_bytes: int = 0
    ui_sent: int = 0
    starts_sent: int = 0
    stops_sent: int = 0

    def note_first_reply(self, ms: float) -> None:
        self.first_reply_lat_ms_accum += ms
        self.first_reply_lat_count += 1

    @property
    def avg_first_reply_ms(self) -> Optional[float]:
        if self.first_reply_lat_count == 0:
            return None
        return self.first_reply_lat_ms_accum / self.first_reply_lat_count


# ---- error capture -----------------------------------------------------------


@dataclass
class ErrorEvent:
    session_id: int
    phase: (
        str  # 'handshake' | 'send' | 'recv' | 'parse' | 'close' | 'no_reply' | 'other'
    )
    when_s: float  # seconds since program start
    exc_type: str
    message: str
    detail: Optional[str] = None  # traceback if --trace enabled


# ---- helpers: build protobuf messages ---------------------------------------


def _mk_user_interaction(seq, Z, xyz, vel=None, cell_size=15.0):
    m = pb.ClientAction()
    m.schema_version = 1
    m.seq = int(seq)
    ui = pb.ClientAction.UserInteraction()
    ui.atomic_numbers.extend(int(z) for z in Z)
    for p in xyz:
        ui.positions.extend([float(p[0]), float(p[1]), float(p[2])])
    if vel is not None:
        for v in vel:
            ui.velocities.extend([float(v[0]), float(v[1]), float(v[2])])
    cell = pb.Mat3()
    cell.m.extend([cell_size, 0.0, 0.0, 0.0, cell_size, 0.0, 0.0, 0.0, cell_size])
    ui.cell.CopyFrom(cell)
    m.user_interaction.CopyFrom(ui)
    return m


def _mk_start(seq, mode, *, calculator):
    m = pb.ClientAction()
    m.schema_version = 1
    m.seq = int(seq)
    start = pb.ClientAction.Start()
    start.simulation_type = (
        pb.ClientAction.Start.SimType.MD
        if str(mode).lower() == "md"
        else pb.ClientAction.Start.SimType.RELAX
    )
    sp = start.simulation_params
    sp.calculator = str(calculator)
    sp.temperature = 298.0
    sp.timestep_fs = 1.0
    sp.friction = 0.02
    sp.fmax = 0.05
    sp.max_step = 0.2
    sp.optimizer = "bfgs"
    m.start.CopyFrom(start)
    return m


def _mk_stop(seq):
    m = pb.ClientAction()
    m.schema_version = 1
    m.seq = int(seq)
    m.stop.CopyFrom(pb.ClientAction.Stop())
    return m


def _mk_ping(seq, ack):
    m = pb.ClientAction()
    m.schema_version = 1
    m.seq = int(seq)
    m.ack = int(ack)
    m.ping.CopyFrom(pb.ClientAction.Ping())
    return m


# ---- debug helpers -----------------------------------------------------------


def _peek(b: bytes, n: int) -> str:
    try:
        return b[:n].decode("utf-8", errors="replace")
    except Exception:
        return str(b[:n])


def _sniff_json_bytes(b: bytes):
    try:
        i = 0
        while i < len(b) and b[i] <= 32:
            i += 1
        if i < len(b) and b[i] in (ord("{"), ord("[")):
            return json.loads(b.decode("utf-8", errors="ignore"))
    except Exception:
        pass
    return None


def _jitter_xyz(xyz, sigma=0.02):
    import random as _r

    return [
        [
            p[0] + _r.gauss(0, sigma),
            p[1] + _r.gauss(0, sigma),
            p[2] + _r.gauss(0, sigma),
        ]
        for p in xyz
    ]


def _rand_vel(n, scale=0.01):
    import random as _r

    return [
        [
            _r.uniform(-scale, scale),
            _r.uniform(-scale, scale),
            _r.uniform(-scale, scale),
        ]
        for _ in range(n)
    ]


# ---- session task -----------------------------------------------------------


async def _one_session(
    session_id: int,
    start_ts: float,
    uri: str,
    m: Metrics,
    stop_time: float,
    debug: bool,
    peek_bytes: int,
    # mix controls
    ui_rate: float,
    start_prob: float,
    stop_prob: float,
    min_run_sec: float,
    md_prob: float,
    calculator: str,
    # error sink
    err_list: List[ErrorEvent],
    include_trace: bool,
    open_timeout: float,
):
    def now_s() -> float:
        return time.perf_counter() - start_ts

    def dlog(*a):
        if debug:
            print(f"[S{session_id}]", *a, flush=True)

    def note_error(
        phase: str, e_type: str, msg: str, exc: Optional[BaseException] = None
    ):
        evt = ErrorEvent(
            session_id=session_id,
            phase=phase,
            when_s=now_s(),
            exc_type=e_type,
            message=msg,
            detail=(tb.format_exc() if (include_trace and exc is not None) else None),
        )
        err_list.append(evt)
        print(
            f"[S{session_id}] ERROR[{phase}] {evt.exc_type}: {evt.message}", flush=True
        )
        if include_trace and evt.detail and debug:
            print(evt.detail, flush=True)

    # Initial seed
    Z, xyz0 = _water()
    xyz = [p[:] for p in xyz0]

    first_reply_at = None
    t_start = time.perf_counter()
    last_seq = 0
    seq = 1
    ack_period = 0.10

    running = False
    current_mode = None  # "md" or "relax"
    last_start_t = 0.0

    def schedule_next_ui(now):
        if ui_rate <= 0:
            return now + 9999.0
        import random as _r

        return now + _r.expovariate(ui_rate)

    try:
        # CONNECT/HANDSHAKE (explicit open_timeout)
        try:
            ws = await websockets.connect(
                uri,
                max_queue=None,
                compression=None,
                ping_interval=20,
                ping_timeout=20,
                max_size=None,
                open_timeout=open_timeout,
            )
        except Exception as e:
            note_error("handshake", type(e).__name__, str(e), e)
            m.sessions_err += 1
            return

        try:
            async with ws:
                dlog("connected", uri)

                # Prime the server with an initial UI
                try:
                    vel = _rand_vel(len(Z), 0.0)  # zero velocities for first push
                    ui = _mk_user_interaction(seq, Z, xyz, vel)
                    await ws.send(ui.SerializeToString())
                    dlog("tx UI(seq=", seq, ") natoms=", len(Z))
                    m.ui_sent += 1
                    seq += 1
                except Exception as e:
                    note_error("send", type(e).__name__, str(e), e)
                    m.sessions_err += 1
                    return

                next_ui_at = schedule_next_ui(time.perf_counter())

                # Acknowledge proactively
                try:
                    ping0 = _mk_ping(seq, 0)
                    await ws.send(ping0.SerializeToString())
                    dlog("tx ack (proactive) seq=", seq, "ack=0")
                    seq += 1
                except Exception as e:
                    note_error("send", type(e).__name__, str(e), e)
                    m.sessions_err += 1
                    return

                last_ack_sent_at = 0.0

                while time.perf_counter() < stop_time:
                    timeout = max(0.01, min(0.10, stop_time - time.perf_counter()))
                    if timeout <= 0:
                        break

                    now = time.perf_counter()

                    # USER_INTERACTION timer
                    if now >= next_ui_at:
                        try:
                            xyz = _jitter_xyz(xyz, sigma=0.02)
                            include_vel = random.random() < 0.5
                            vel = _rand_vel(len(Z), 0.02) if include_vel else None
                            ui = _mk_user_interaction(seq, Z, xyz, vel)
                            await ws.send(ui.SerializeToString())
                            dlog("tx UI(seq=", seq, ") vel=", bool(include_vel))
                            m.ui_sent += 1
                            seq += 1
                            next_ui_at = schedule_next_ui(now)
                        except Exception as e:
                            note_error("send", type(e).__name__, str(e), e)
                            break

                    # Start
                    if (not running) and random.random() < start_prob:
                        try:
                            current_mode = (
                                "md" if (random.random() < md_prob) else "relax"
                            )
                            st = _mk_start(seq, current_mode, calculator=calculator)
                            await ws.send(st.SerializeToString())
                            dlog("tx START seq=", seq, "mode=", current_mode)
                            m.starts_sent += 1
                            seq += 1
                            running = True
                            last_start_t = now
                            ping = _mk_ping(seq, last_seq)
                            await ws.send(ping.SerializeToString())
                            dlog("tx ack (post-start) seq=", seq, "ack=", last_seq)
                            seq += 1
                        except Exception as e:
                            note_error("send", type(e).__name__, str(e), e)
                            break

                    # Stop (only after min_run_sec)
                    if (
                        running
                        and (now - last_start_t >= min_run_sec)
                        and (random.random() < stop_prob)
                    ):
                        try:
                            stop = _mk_stop(seq)
                            await ws.send(stop.SerializeToString())
                            dlog("tx STOP seq=", seq)
                            m.stops_sent += 1
                            seq += 1
                            running = False
                            current_mode = None
                        except Exception as e:
                            note_error("send", type(e).__name__, str(e), e)
                            break

                    # Receive
                    try:
                        msg = await asyncio.wait_for(ws.recv(), timeout=timeout)
                    except asyncio.TimeoutError:
                        now2 = time.perf_counter()
                        if now2 - last_ack_sent_at >= ack_period and last_seq > 0:
                            try:
                                ping = _mk_ping(seq, last_seq)
                                await ws.send(ping.SerializeToString())
                                dlog("tx ack (periodic) seq=", seq, "ack=", last_seq)
                                seq += 1
                                last_ack_sent_at = now2
                            except Exception as e:
                                note_error("send", type(e).__name__, str(e), e)
                                break
                        continue
                    except wse.ConnectionClosedOK as e:
                        # Clean close (OK); treat as ok if we had a reply, else synthesize no_reply
                        if first_reply_at is None:
                            note_error(
                                "no_reply",
                                type(e).__name__,
                                f"closed before first reply (code={getattr(e,'code',None)} reason={getattr(e,'reason',None)})",
                                e,
                            )
                        else:
                            # Some servers close with code 1000 mid-churn; treat as ok
                            dlog("closed OK")
                        break
                    except wse.ConnectionClosedError as e:
                        # Abnormal close
                        note_error(
                            "close",
                            type(e).__name__,
                            f"code={getattr(e,'code',None)} reason={getattr(e,'reason',None)}",
                            e,
                        )
                        break
                    except Exception as e:
                        note_error("recv", type(e).__name__, str(e), e)
                        break

                    # Process message
                    if isinstance(msg, str):
                        m.frames_text += 1
                        if first_reply_at is None:
                            first_reply_at = time.perf_counter()
                        m.messages_rcvd += 1
                        continue

                    if isinstance(msg, (bytes, bytearray, memoryview)):
                        raw = bytes(msg)
                        sr = pb.ServerResult()
                        try:
                            sr.ParseFromString(raw)  # type: ignore
                            m.frames_proto_ok += 1
                        except Exception as e:
                            m.frames_proto_parse_err += 1
                            note_error(
                                "parse",
                                type(e).__name__,
                                f"proto parse failed; peek={_peek(raw, min(len(raw), peek_bytes))}",
                                e,
                            )
                            continue

                        if sr.seq:
                            last_seq = max(last_seq, int(sr.seq))
                        if first_reply_at is None:
                            first_reply_at = time.perf_counter()
                        m.messages_rcvd += 1
                        continue

                # session end
                if first_reply_at is not None:
                    m.sessions_ok += 1
                    m.note_first_reply((first_reply_at - t_start) * 1000.0)
                else:
                    m.sessions_err += 1
                    note_error(
                        "no_reply",
                        "NoFirstReply",
                        "session ended without any server reply",
                        None,
                    )

        except Exception as e:
            # async with body raised
            m.sessions_err += 1
            note_error("other", type(e).__name__, str(e), e)

    except Exception as e:
        m.sessions_err += 1
        note_error("other", type(e).__name__, str(e), e)


# ---- driver -----------------------------------------------------------------


async def run_stress_ws(
    base_ws: str,
    duration_s: int,
    concurrency: int,
    session_mode: str,
    debug: bool,
    peek_bytes: int,
    ui_rate: float,
    start_prob: float,
    stop_prob: float,
    min_run_sec: float,
    md_prob: float,
    calculator: str,
    error_log: Optional[str],
    include_trace: bool,
    open_timeout: float,
):
    m = Metrics()
    base = base_ws.rstrip("/") + "/ws"
    started = time.perf_counter()
    deadline = started + float(duration_s)

    if debug:
        print(
            f"[driver] base={base} duration={duration_s}s concurrency={concurrency} session_mode={session_mode} ui_rate={ui_rate} start_prob={start_prob} stop_prob={stop_prob} min_run_sec={min_run_sec} md_prob={md_prob} open_timeout={open_timeout}",
            flush=True,
        )

    tasks = []
    sid = 1
    err_list: List[ErrorEvent] = []

    async def ticker_task():
        i = 0
        try:
            while time.perf_counter() < deadline:
                await asyncio.sleep(1.0)
                i += 1
                print(
                    f"[tick {i:02d}s] msgs={m.messages_rcvd} ok={m.sessions_ok} err={m.sessions_err} proto_ok={m.frames_proto_ok} proto_err={m.frames_proto_parse_err} text={m.frames_text} json_bytes={m.frames_json_in_bytes} other_bytes={m.frames_other_bytes} ui_sent={m.ui_sent} starts={m.starts_sent} stops={m.stops_sent}",
                    flush=True,
                )
        except asyncio.CancelledError:
            pass

    t_task = asyncio.create_task(ticker_task())

    try:
        if session_mode == "steady":
            stop_time = deadline
            for _ in range(concurrency):
                tasks.append(
                    asyncio.create_task(
                        _one_session(
                            sid,
                            started,
                            base,
                            m,
                            stop_time,
                            debug,
                            peek_bytes,
                            ui_rate,
                            start_prob,
                            stop_prob,
                            min_run_sec,
                            md_prob,
                            calculator,
                            err_list,
                            include_trace,
                            open_timeout,
                        )
                    )
                )
                sid += 1
            await asyncio.gather(*tasks, return_exceptions=True)
        else:
            while time.perf_counter() < deadline:
                tasks = [t for t in tasks if not t.done()]
                while len(tasks) < concurrency and time.perf_counter() < deadline:
                    stop_time = min(deadline, time.perf_counter() + (min_run_sec + 2.0))
                    tasks.append(
                        asyncio.create_task(
                            _one_session(
                                sid,
                                started,
                                base,
                                m,
                                stop_time,
                                debug,
                                peek_bytes,
                                ui_rate,
                                start_prob,
                                stop_prob,
                                min_run_sec,
                                md_prob,
                                calculator,
                                err_list,
                                include_trace,
                                open_timeout,
                            )
                        )
                    )
                    sid += 1
                await asyncio.sleep(0.01)
            if tasks:
                await asyncio.gather(*tasks, return_exceptions=True)
    finally:
        t_task.cancel()
        try:
            await t_task
        except Exception:
            pass

    # Optional JSONL error log
    if error_log and err_list:
        try:
            with open(error_log, "w", encoding="utf-8") as f:
                for e in err_list:
                    f.write(json.dumps(asdict(e), ensure_ascii=False) + "\n")
            print(f"[errors] wrote {len(err_list)} errors to {error_log}", flush=True)
        except Exception as e:
            print(f"[errors] failed to write error log: {e}", flush=True)

    return m, err_list


# ---- error summary printing --------------------------------------------------


def _print_error_summary(errs: List[ErrorEvent], include_trace: bool):
    if not errs:
        print("\nNo errors captured.")
        return
    print("\n--- Error Summary ---")
    print(f"Total errors: {len(errs)}")
    # by phase
    by_phase: Dict[str, int] = {}
    by_exc: Dict[str, int] = {}
    for e in errs:
        by_phase[e.phase] = by_phase.get(e.phase, 0) + 1
        key = f"{e.exc_type}"
        by_exc[key] = by_exc.get(key, 0) + 1
    print("By phase: " + ", ".join(f"{k}={v}" for k, v in sorted(by_phase.items())))
    print("By exception: " + ", ".join(f"{k}={v}" for k, v in sorted(by_exc.items())))
    # first few detailed examples
    print("\nExamples:")
    for e in errs[: min(10, len(errs))]:
        print(
            f"  [S{e.session_id}] t=+{e.when_s:.2f}s phase={e.phase} {e.exc_type}: {e.message}"
        )
        if include_trace and e.detail:
            print("  Traceback:")
            print("  " + "\n  ".join(line.rstrip() for line in e.detail.splitlines()))


# ---- CLI --------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="UMA WS stress utility (MIXED+ERRORS plus)"
    )
    ap.add_argument("--base", default="ws://127.0.0.1:8000")
    ap.add_argument("--duration", type=int, default=30)
    ap.add_argument("--concurrency", type=int, default=8)
    ap.add_argument("--session-mode", choices=["steady", "churn"], default="churn")
    ap.add_argument("--debug", action="store_true")
    ap.add_argument("--peek-bytes", type=int, default=96)
    ap.add_argument("--ui-rate", type=float, default=8.0)
    ap.add_argument("--start-prob", type=float, default=0.12)
    ap.add_argument("--stop-prob", type=float, default=0.06)
    ap.add_argument("--min-run-sec", type=float, default=2.0)
    ap.add_argument("--md-prob", type=float, default=0.7)
    ap.add_argument("--calculator", choices=["uma", "lj"], default="uma")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--error-log", type=str, default=None)
    ap.add_argument("--trace", action="store_true")
    ap.add_argument(
        "--open-timeout",
        type=float,
        default=20.0,
        help="WS opening handshake timeout (seconds)",
    )

    args = ap.parse_args()
    if args.seed is not None:
        random.seed(args.seed)

    m, errs = asyncio.run(
        run_stress_ws(
            base_ws=args.base,
            duration_s=args.duration,
            concurrency=args.concurrency,
            session_mode=args.session_mode,
            debug=args.debug,
            peek_bytes=args.peek_bytes,
            ui_rate=args.ui_rate,
            start_prob=args.start_prob,
            stop_prob=args.stop_prob,
            min_run_sec=args.min_run_sec,
            md_prob=args.md_prob,
            calculator=args.calculator,
            error_log=args.error_log,
            include_trace=args.trace,
            open_timeout=args.open_timeout,
        )
    )

    print("\n=== UMA WS Stress Results (MIXED+ERRORS+) ===")
    print(f"base= {args.base}")
    print(
        f"duration_s= {args.duration} concurrency= {args.concurrency} session_mode= {args.session_mode}"
    )
    print(f"sessions_ok= {m.sessions_ok} sessions_err= {m.sessions_err}")
    print(f"messages_rcvd= {m.messages_rcvd}")
    print(
        f"frames_proto_ok= {m.frames_proto_ok} frames_proto_parse_err= {m.frames_proto_parse_err}"
    )
    print(
        f"frames_text= {m.frames_text} frames_json_in_bytes= {m.frames_json_in_bytes} frames_other_bytes= {m.frames_other_bytes}"
    )
    print(
        f"ui_sent= {m.ui_sent} starts_sent= {m.starts_sent} stops_sent= {m.stops_sent}"
    )
    if m.avg_first_reply_ms is not None:
        print(f"avg_first_reply_ms= {m.avg_first_reply_ms:.2f}")
    else:
        print("avg_first_reply_ms= n/a")

    _print_error_summary(errs, include_trace=args.trace)


if __name__ == "__main__":
    main()
