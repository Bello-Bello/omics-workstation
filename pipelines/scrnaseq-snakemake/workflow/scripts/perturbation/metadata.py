"""Perturbation condition parsing shared by the Nextflow/Snakemake/notebook entry points.

Kept separate from qc.py/analysis.py because condition-column naming is the
one thing that genuinely differs between perturbation datasets (sciPlex uses
"perturbation"/"dose_value", other screens use "condition"/"dose", etc.) —
isolating the assumption here means qc.py and analysis.py can stay dataset-agnostic.
"""

from __future__ import annotations

import anndata as ad
import pandas as pd


def standardize_condition_columns(
    adata: ad.AnnData,
    perturbation_col: str,
    dose_col: str | None,
    control_label: str,
) -> ad.AnnData:
    """Rename dataset-specific condition columns to a fixed schema.

    Downstream steps key off `adata.obs["condition"]` and
    `adata.obs["is_control"]` regardless of what the source dataset called
    them, so this is the only place that needs to change to point the
    pipeline at a different perturbation dataset.
    """
    adata.obs["condition"] = adata.obs[perturbation_col].astype(str)
    adata.obs["is_control"] = adata.obs["condition"] == control_label

    if dose_col is not None:
        # sciPlex encodes dose as a string (categories like "0.1", "100") with
        # NaN for the control wells — coerce to float so dose-response
        # ordering/plotting works, leaving control rows as NaN rather than 0
        # (0 would misleadingly imply "dose zero of a real compound").
        adata.obs["dose"] = pd.to_numeric(adata.obs[dose_col], errors="coerce")

    return adata


def validate_conditions(adata: ad.AnnData) -> dict:
    """Check that every non-control condition has a usable control to compare against.

    A perturbation experiment with a condition but no matching control cells
    (e.g. lost entirely to QC filtering) can't support differential
    expression — better to catch that here than have `rank_genes_groups`
    fail cryptically downstream.
    """
    if "condition" not in adata.obs.columns:
        raise ValueError("adata.obs missing 'condition' — run standardize_condition_columns first")

    n_control = int(adata.obs["is_control"].sum())
    counts = adata.obs.loc[~adata.obs["is_control"], "condition"].value_counts()

    return {
        "n_control_cells": n_control,
        "has_control": n_control > 0,
        "conditions": counts.to_dict(),
    }
