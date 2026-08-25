"""Runs inside Snakemake's `script:` directive — `snakemake` object is injected."""

import sys
from pathlib import Path

import pertpy as pt

DATASETS = {
    "sciplex2": pt.data.srivatsan_2020_sciplex2,
}

name = snakemake.params["dataset"]
out_path = Path(snakemake.output[0])

if name not in DATASETS:
    sys.exit(f"Unknown dataset {name!r}. Available: {list(DATASETS)}")

out_path.parent.mkdir(parents=True, exist_ok=True)
adata = DATASETS[name]()
adata.write_h5ad(out_path)
print(f"{name}: {adata.n_obs} cells x {adata.n_vars} genes -> {out_path}")
