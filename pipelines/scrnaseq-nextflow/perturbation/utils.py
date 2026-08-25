"""Logging and I/O helpers shared across the perturbation package."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Module logger with a consistent format for both notebook and Nextflow/Snakemake log output."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(level)
    return logger


def load_adata(path: str | Path) -> ad.AnnData:
    path = Path(path)
    if path.suffix != ".h5ad":
        raise ValueError(f"Expected a .h5ad file, got: {path}")
    return ad.read_h5ad(path)


def save_adata(adata: ad.AnnData, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
