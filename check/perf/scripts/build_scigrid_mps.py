#!/usr/bin/env python3
"""Build a PyPSA-derived energy-model LP in MPS format for HiGHS benchmarking.

Uses PyPSA's bundled SciGRID-DE example network (~585 buses, ~850 lines,
~1400 generators). The 24h example is too small to be a meaningful end-to-end
solver benchmark, so we replicate its snapshot index over a configurable
horizon to scale up the LP without needing extra source data.

Approximate sizing on a 2.5 GHz x86_64 (HiGHS Release):

    horizon  variables   constraints   mps size   solve time
    -------  ---------   -----------   --------   ----------
       24h     ~60 000      ~143 000     ~14 MB         ~5 s
      168h    ~417 000     ~1.0 M       ~99 MB        ~50 s
      336h    ~835 000     ~2.0 M      ~200 MB       ~2 min
      504h     ~1.25 M     ~3.0 M      ~300 MB     ~5-6 min

The 504h variant is the recommended "real workload" for evaluating CPU
optimisations whose effect on simplex inner-loops only shows up at scale.

Usage:
    python3 build_scigrid_mps.py --hours 504 --out /tmp/scigrid_504h.mps
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import pandas as pd
import pypsa


def build(hours: int, out_path: str) -> None:
    t0 = time.time()
    n = pypsa.examples.scigrid_de(from_master=False)
    n.set_snapshots(pd.date_range(start="2011-01-01", periods=hours, freq="h"))
    print(
        f"[scigrid] network: {len(n.buses)} buses, {len(n.lines)} lines, "
        f"{len(n.generators)} generators, {hours} snapshots",
        file=sys.stderr,
    )
    m = n.optimize.create_model()
    print(
        f"[scigrid] linopy model: nvars={m.nvars} ncons={m.ncons}",
        file=sys.stderr,
    )
    m.to_file(out_path)
    size_mb = os.path.getsize(out_path) / 1e6
    print(
        f"[scigrid] wrote {out_path} ({size_mb:.1f} MB) in "
        f"{time.time() - t0:.1f}s",
        file=sys.stderr,
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--hours",
        type=int,
        default=504,
        help="Snapshot horizon (default 504h ~ 3 weeks; ~5-6 min on HiGHS).",
    )
    p.add_argument(
        "--out",
        default="/tmp/scigrid.mps",
        help="Output MPS path (default /tmp/scigrid.mps).",
    )
    args = p.parse_args()
    build(args.hours, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
