import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

try:
    from ase.build import molecule
    from ase.optimize import BFGS
    from ase.calculators.lj import LennardJones
except Exception:  # pragma: no cover
    molecule = None
    BFGS = None
    LennardJones = None


SERVER_URL = os.environ.get("RELAX_PARITY_SERVER", "http://127.0.0.1:8000")
REPO_ROOT = Path(__file__).resolve().parents[1]
NODE = os.environ.get("NODE_BIN", "node")
SCRIPT = REPO_ROOT / "scripts" / "relax_via_api.js"


def have_server():
    """Return True if the FastAPI server responds 200 to a minimal valid request.

    We send a syntactically valid /simple_calculate body (1 atom, dummy coords).
    A 422 previously caused false negatives; only treat network/5xx as failure.
    """
    import urllib.request
    import json as _json
    payload = _json.dumps({
        "atomic_numbers": [1],
        "coordinates": [[0,0,0]],
        "properties": ["energy"],
        "calculator": "lj"
    }).encode()
    req = urllib.request.Request(
        f"{SERVER_URL}/simple_calculate",
        method="POST",
        data=payload,
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=2) as resp:
            code = getattr(resp, 'status', resp.getcode())
            return 200 <= code < 300
    except Exception:
        return False


@pytest.mark.parametrize("calculator", ["lj", "uma"])
def test_js_relax_matches_direct(calculator):
    if not have_server():
        pytest.skip("Server not reachable on %s" % SERVER_URL)
    if molecule is None or BFGS is None:
        pytest.skip("ASE not installed in test environment")

    # Build water (ASE ordering H2O => H,H,O; we reorder to O,H,H for consistency with frontend)
    at = molecule("H2O")
    # Reorder to O,H,H (frontend typical ordering) and translate O to origin
    # ase.build.molecule puts O last. Capture positions.
    pos = at.get_positions()
    # indices: 0 H,1 H,2 O -> reorder to 2,0,1
    reorder = [2, 0, 1]
    pos_re = pos[reorder]
    pos_re[:, :] -= pos_re[0]  # shift O to origin
    # Ensure plain Python ints (avoid numpy.int64 JSON serialization issues)
    numbers = [int(at.numbers[i]) for i in reorder]

    # Direct ASE relaxation (fixed steps = 20 per parity requirement)
    at2 = molecule("H2O")
    at2 = at2[[2, 0, 1]]  # reorder
    at2.positions[:] = pos_re
    if calculator == "lj":
        if LennardJones is None:
            pytest.skip("ASE LennardJones unavailable")
        at2.calc = LennardJones(rc=3.0)
    else:  # uma -> exercise /simple_calculate path would require FairChem; skip if unavailable
        try:
            from fairchem.core import pretrained_mlip
            # Updated path for FAIRChemCalculator in current fairchem version
            from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
            from fairchem.core.units.mlip_unit import InferenceSettings
        except Exception:
            pytest.skip("FairChem not available for UMA direct reference")
        pu = pretrained_mlip.get_predict_unit(
            os.environ.get("UMA_MODEL", "uma-s-1p1"),
            device="cpu",
            inference_settings=InferenceSettings(
                tf32=False,
                activation_checkpointing=False,
                merge_mole=True,
                compile=False,
                external_graph_gen=False,
                internal_graph_gen_version=2,
            ),
        )
        at2.calc = FAIRChemCalculator(pu, task_name=os.environ.get("UMA_TASK", "omol"))

    initial_energy_ref = float(at2.get_potential_energy())
    opt = BFGS(at2, logfile=None)
    steps = 20
    for _ in range(steps):
        opt.step()
    final_energy_ref = float(at2.get_potential_energy())

    # Call server via Node script (JS client path)
    cmd = [
        NODE,
        str(SCRIPT),
        "--url",
        SERVER_URL,
        "--calc",
        calculator,
        "--steps",
        str(steps),
        "--atomic-numbers",
        json.dumps(numbers),
        "--coords",
        json.dumps(pos_re.tolist()),
    ]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    assert proc.returncode == 0, f"Node script failed: {proc.stderr}"
    data = json.loads(proc.stdout)
    assert data["steps_completed"] == steps

    # Energy comparisons: allow slightly looser tolerance for UMA due to ML noise vs different graph init
    if calculator == "lj":
        assert abs(data["initial_energy"] - initial_energy_ref) < 1e-4
        # Lennard-Jones deterministic; final should match tightly
        assert abs(data["final_energy"] - final_energy_ref) < 1e-4
    else:
        # UMA: tolerate small drift (model caching differences, device, etc.)
        assert abs(data["initial_energy"] - initial_energy_ref) < 5e-3
        assert abs(data["final_energy"] - final_energy_ref) < 5e-3

    # Shape checks
    assert len(data["positions"]) == len(numbers)
    assert len(data["forces"]) == len(numbers)
    assert all(len(f) == 3 for f in data["forces"])  # 3D forces
    print(
        f"[parity] calc={calculator} dt={dt*1000:.1f}ms E0_ref={initial_energy_ref:.6f} E0_api={data['initial_energy']:.6f} Ef_ref={final_energy_ref:.6f} Ef_api={data['final_energy']:.6f}"
    )
