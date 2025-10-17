from __future__ import annotations

import argparse
import time

from fairchem_local_server2.ws_app import deploy


def main():
    parser = argparse.ArgumentParser(description="Start UMA WS Serve app")
    parser.add_argument("--ngpus", type=int, default=None)
    parser.add_argument("--ncpus", type=int, default=None)
    parser.add_argument("--nhttp", type=int, default=None)
    args = parser.parse_args()

    deploy(args.ngpus, args.ncpus, args.nhttp)
    try:
        while True:
            time.sleep(60)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
