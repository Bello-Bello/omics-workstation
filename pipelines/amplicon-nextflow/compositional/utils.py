"""Shared helpers for the compositional analysis package."""

from __future__ import annotations

import logging
import sys


def get_logger(name: str) -> logging.Logger:
    """Logger that writes to stderr so it never contaminates piped stdout.

    Matches the convention used by the perturbation package in
    ../scrnaseq-nextflow/perturbation/utils.py — Nextflow captures both
    streams, but keeping stdout clean means a process can still be swapped
    to a `... > out.tsv` form without the log ending up in the data.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger
