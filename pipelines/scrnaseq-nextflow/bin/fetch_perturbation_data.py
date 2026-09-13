#!/usr/bin/env python3
"""Fetch a public perturbation screen dataset and cache it as .h5ad.

Called as: fetch_perturbation_data.py <dataset_name> <output_path>

Only one dataset is wired up (sciplex2) rather than a generic loader —
adding a second real dataset later is a one-line addition to DATASETS, which
does not justify a plugin system at this scale.
"""

import sys
from pathlib import Path

import pertpy as pt

DATASETS = {
    "sciplex2": pt.data.srivatsan_2020_sciplex2,
}


def main() -> None:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <dataset_name> <output_path>", file=sys.stderr)
        sys.exit(1)

    name, out_path = sys.argv[1], Path(sys.argv[2])
    if name not in DATASETS:
        print(f"Unknown dataset {name!r}. Available: {list(DATASETS)}", file=sys.stderr)
        sys.exit(1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata = DATASETS[name]()
    adata.write_h5ad(out_path)
    print(f"{name}: {adata.n_obs} cells x {adata.n_vars} genes -> {out_path}")


if __name__ == "__main__":
    main()
