#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute TASK-CONSTRAINED INTRA-HEMISPHERIC PARTIAL correlations between
ROI macro-regions (seeds) and TASK-DEFINED macro-region targets using
Nilearn ConnectivityMeasure.

Differs from the parcel-level (task-free) script in two key respects:

1. Targets are task-defined macro-regions (same binary masks as seeds,
   loaded from the /target/ folder), not atlas MMP parcels.
   Output matrix shape: (n_clusters × n_clusters).

2. Conditioning set = the OTHER macro-region timeseries within the same
   task-defined network (excluding seed and target).
   This mirrors the logic of the task-free script exactly:
     - Task-free:        condition on all MMP parcels except seed's own + target
     - Task-constrained: condition on all task macro-regions except seed + target
   Full vs partial correlations are therefore directly comparable within
   each analysis variant (task-free / task-constrained).

   Note: conditioning on atlas parcels outside the task network would change
   what is being measured and break the symmetry with the task-free analysis.

Design decisions
----------------
Seed timeseries   : ipsilateral hemisphere only (macro-region mask → mean vertex signal)
Target timeseries : ipsilateral (primary output) + contralateral (bilateral output)
Conditioning set  : remaining task macro-regions, ipsilateral only
                    (same hemisphere restriction as in the task-free script)

Self-correlation masking
  Seed == target entries are left as NaN (diagonal of the output matrix).

Standardize flag
  'zscore_sample' throughout: safe on already-z-scored XCP-D data (near-no-op
  since mean≈0, std≈1) and correct on non-z-scored individual runs.
  Forward-compatible with Nilearn ≥ 0.15.

Collinearity diagnostic + regularized covariance estimator (UPDATED)
  A first diagnostic pass (condition number of the covariance matrix of
  the full conditioning-relevant region pool, computed per subject × run
  × hemisphere) showed condition numbers in the hundreds to thousands for
  every single subject/run/hemi in this dataset (120/120 flagged HIGH,
  median ~1400, max ~8400).  This indicates the raw sample covariance
  matrix is close to singular, which is the classic driver of unstable /
  high-variance partial correlation estimates (inverting a near-singular
  matrix amplifies noise disproportionately).

  Given this, ConnectivityMeasure now uses cov_estimator=LedoitWolf()
  instead of the default (unregularized) empirical covariance.  Ledoit-
  Wolf shrinkage (Ledoit & Wolf, 2004) shrinks the sample covariance
  matrix toward a scaled identity matrix by an amount computed
  analytically from the data — no cross-validation or hyperparameter
  tuning required — which stabilizes inversion without changing which
  regions are in the conditioning set (so full-corr vs. partial-corr
  comparability is preserved).  Conceptually related to the "graphical
  ridge" estimator benchmarked in Peterson et al. (2025, Imaging
  Neuroscience), which found regularized partial correlation improves
  reliability over the unregularized estimate in fMRI FC.

  The condition-number diagnostic is still computed for every
  subject/run/hemi (on the RAW covariance, as before) so the shrinkage
  benefit can be inspected quantitatively: after computing the raw
  condition number, the same check is repeated on the covariance matrix
  AFTER Ledoit-Wolf shrinkage, and both are logged side by side.  Neither
  diagnostic changes the saved partial-corr outputs — only cov_estimator
  does.

Outputs (per subject, per run, per hemisphere)
  seed-task_by_macror-task_partial_{run}_{hemi}.npy / .csv           — Pearson r
  seed-task_by_macror-task_partial_fisherz_{run}_{hemi}.npy / .csv   — Fisher z
  seed-task_by_macror-task_partial_{run}_{hemi}_bilateral.npy / .csv — bilateral

Diagnostic output (all subjects/runs/hemis combined)
  {main_data}/group/91k/rest/partial_corr/diagnostics/
      collinearity_condition_numbers_task-constrained.csv

Group aggregation: run group_stats_partial_corr_task-constrained.py 
(Fisher-z space throughout)

Reference: Peterson, K.L., Sanchez-Romero, R., Mill, R.D., & Cole, M.W. (2025).
Regularized partial correlation provides reliable functional connectivity
estimates while correcting for widespread confounding. Imaging Neuroscience.
https://doi.org/10.1162/IMAG.a.162

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd
from nilearn.connectome import ConnectivityMeasure
from sklearn.covariance import LedoitWolf, EmpiricalCovariance, GraphicalLassoCV

# ============================================================
# Covariance estimator
#
# ESTIMATOR_TAG is passed in as sys.argv[1] by the SLURM submit script
# (submit_nilearn_compute_partial_corr_job.py), which owns the choice of
# estimator as a run-configuration decision. This script only validates
# the tag and maps it to the corresponding sklearn estimator — it does
# NOT default silently, since running without an explicit choice on the
# cluster would make it too easy to lose track of which variant a given
# job produced.
#
# raw            : EmpiricalCovariance (unregularized) — Nilearn's default
# ledoit-wolf    : LedoitWolf shrinkage (Ledoit & Wolf, 2004) — analytic
#                  shrinkage intensity, no cross-validation
# graphical-lasso: GraphicalLassoCV — L1-regularized precision matrix,
#                  sparsity penalty selected by cross-validation. Slower
#                  than the other two; produces a sparse (many-zero)
#                  precision matrix, which is a different scientific claim
#                  than shrinkage (see Peterson et al. 2025, Imaging
#                  Neuroscience, for the graphical lasso vs. graphical
#                  ridge distinction).
#
# Output filenames are tagged with ESTIMATOR_TAG so that runs with
# different estimators never overwrite each other on disk.
# ============================================================

VALID_ESTIMATORS = ("raw", "ledoit-wolf", "graphical-lasso")

if len(sys.argv) < 2:
    print(
        "ERROR: no covariance estimator specified.\n"
        f"  Usage: python {sys.argv[0]} <estimator>\n"
        f"  Accepted: {', '.join(VALID_ESTIMATORS)}\n"
        "  This should normally be set via submit_nilearn_compute_partial_corr_job.py, "
        "not called directly."
    )
    sys.exit(1)

ESTIMATOR_TAG = sys.argv[1]
if ESTIMATOR_TAG not in VALID_ESTIMATORS:
    print(
        f"ERROR: unrecognised estimator '{ESTIMATOR_TAG}'.\n"
        f"  Accepted: {', '.join(VALID_ESTIMATORS)}"
    )
    sys.exit(1)

if ESTIMATOR_TAG == "ledoit-wolf":
    COV_ESTIMATOR = LedoitWolf()
elif ESTIMATOR_TAG == "graphical-lasso":
    COV_ESTIMATOR = GraphicalLassoCV()
else:
    COV_ESTIMATOR = EmpiricalCovariance()

print(f"Covariance estimator: {ESTIMATOR_TAG}  (cov_estimator={COV_ESTIMATOR})")

# ============================================================
# Paths
# ============================================================

main_data    = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))

sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import impute_nan_columns

# ============================================================
# Settings
# ============================================================

project_dir       = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# ROIs
# ============================================================

clusters = analysis_info["rois-drawn"]
seed_to_parcels = analysis_info["rois-group-mmp"]

# Reverse so mPCS is first (matches downstream visualisation scripts)
clusters.reverse()

# ============================================================
# Hemisphere configuration
# ============================================================

HEMIS = [
    {"label": "LH", "seed_key": "lh", "atlas_key": "L", "ts_key": "data_L"},
    {"label": "RH", "seed_key": "rh", "atlas_key": "R", "ts_key": "data_R"},
]

# ============================================================
# Load rest-specific settings
# ============================================================

rest_settings_path = os.path.join(base_dir, project_dir, "rest-settings.yml")
rest_settings      = load_settings([rest_settings_path])[0]
RUNS          = rest_settings["runs"]

# ============================================================
# Collinearity diagnostic — accumulator
#
# One record per subject × run × hemisphere, capturing the condition
# number of the covariance matrix of the full region pool used for
# conditioning (ipsi + contra macro-regions stacked together), plus
# the number of timepoints and regions for context.  Printed as a
# summary table and saved to CSV at the very end of the script.
# ============================================================
diagnostic_records = []


def _flag_cond(cond):
    if cond < 30:
        return "OK"
    elif cond < 100:
        return "WATCH"
    else:
        return "HIGH"


def _condition_number_report(X, subject, run_tag, label, n_time):
    """
    Compute and log the condition number of the covariance matrix of X
    (n_time × n_regions), both RAW (empirical covariance, same estimator
    Nilearn would use by default) and AFTER Ledoit-Wolf shrinkage (the
    estimator actually used for the partial-corr computation below).
    Appends a record to diagnostic_records but does not alter any
    partial-corr computation — purely informational.
    """
    n_regions = X.shape[1]

    cov_raw   = np.cov(X, rowvar=False)
    cond_raw  = np.linalg.cond(cov_raw)

    cov_lw    = LedoitWolf().fit(X).covariance_
    cond_lw   = np.linalg.cond(cov_lw)

    flag_raw = _flag_cond(cond_raw)
    flag_lw  = _flag_cond(cond_lw)

    print(
        f"  [{label}] Collinearity check: {n_regions} regions, "
        f"{n_time} timepoints — "
        f"RAW cond(cov) = {cond_raw:.1f} [{flag_raw}]  →  "
        f"Ledoit-Wolf cond(cov) = {cond_lw:.2f} [{flag_lw}]"
    )

    diagnostic_records.append({
        "subject":       subject,
        "run":           run_tag,
        "hemi":          label,
        "n_regions":     n_regions,
        "n_timepoints":  n_time,
        "condition_number_raw":         cond_raw,
        "flag_raw":                     flag_raw,
        "condition_number_ledoit_wolf": cond_lw,
        "flag_ledoit_wolf":             flag_lw,
    })


# ============================================================
# Subject loop
# ============================================================

for subject in subjects:

    for run in RUNS:

        run_tag      = f"_{run}" if run else "_concat"   # "" → "_concat" in filenames
        run_tag_ts   = f"_{run}" if run else ""          # "" for timeseries filename (no tag)

        print(f"\n=== Processing {subject}{run_tag} ===")
        
        timeseries_fn = (
            f"{main_data}/{subject}/91k/rest/timeseries/"
            f"{subject}_ses-01_task-rest{run_tag_ts}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
        )

        ts_img, ts_data_raw = load_surface(timeseries_fn)

        res = from_91k_to_32k(
            ts_img, ts_data_raw,
            return_concat_hemis=False,
            return_32k_mask=True,
        )

        mask_32k = res["mask_32k"]
        print(f"  mask_32k: {int(np.sum(~mask_32k))} medial-wall vertices per hemi")

        # ----------------------------------------------------------
        # Hemisphere loop
        # ----------------------------------------------------------

        for h in HEMIS:
            label      = h["label"]     # "LH" or "RH"
            seed_key   = h["seed_key"]  # "lh" or "rh"
            atlas_key  = h["atlas_key"] # "L"  or "R"
            ts_key     = h["ts_key"]    # "data_L" or "data_R"
            contra_key = "rh" if seed_key == "lh" else "lh"
            ts_key_contra = "data_R" if ts_key == "data_L" else "data_L"

            ts_data       = res[ts_key]        # (n_time, n_vertices_ipsi)
            ts_data_contra = res[ts_key_contra] # (n_time, n_vertices_contra)
            print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

            # ------------------------------------------------------
            # Step 1 — Load ALL task macro-region timeseries
            #
            # Seeds, ipsilateral targets, and contralateral targets are all
            # loaded up front into a single dict so the conditioning set
            # can be assembled cleanly in Step 4.
            #
            # Both /seed/ and /target/ use the same file naming convention;
            # they are loaded from separate folders to keep the analysis
            # intent explicit.
            # ------------------------------------------------------

            # macro_ts[roi] = ipsilateral timeseries for that macro-region
            macro_ts        = {}
            macro_ts_contra = {}   # contralateral timeseries
            loaded_rois        = []
            loaded_rois_contra = []

            for roi in clusters:

                # Ipsilateral (used as seed and ipsilateral conditioning)
                _, mask_data = load_surface(
                    f"{main_data}/{subject}/91k/rest/seed/"
                    f"{subject}_91k_intertask_Sac-Pur-pRF_{seed_key}_{roi}.shape.gii"
                )
                mask = mask_data.ravel()
                if np.any(mask):
                    macro_ts[roi] = ts_data[:, mask > 0].mean(axis=1)
                    loaded_rois.append(roi)
                else:
                    print(f"  [{label}] ⚠️  Empty ipsi mask for {roi} — skipping")

                # Contralateral (used as contra target and contralateral conditioning)
                _, mask_data = load_surface(
                    f"{main_data}/{subject}/91k/rest/target/"
                    f"{subject}_91k_intertask_Sac-Pur-pRF_{contra_key}_{roi}.shape.gii"
                )
                mask = mask_data.ravel()
                if np.any(mask):
                    macro_ts_contra[roi] = ts_data_contra[:, mask > 0].mean(axis=1)
                    loaded_rois_contra.append(roi)
                else:
                    print(f"  [{label}] ⚠️  Empty contra mask for {roi} — skipping")

            if not loaded_rois:
                print(f"  [{label}] ⚠️  No valid macro-region timeseries — skipping hemisphere")
                continue

            print(f"  [{label}] Ipsi macro-regions:  {loaded_rois}")
            print(f"  [{label}] Contra macro-regions: {loaded_rois_contra}")

            # ------------------------------------------------------
            # Step 2 — NaN imputation
            # ------------------------------------------------------

            for roi in list(macro_ts.keys()):
                clean = impute_nan_columns(
                    macro_ts[roi][:, np.newaxis],
                    label=f"{subject}{run_tag} {label} {roi}"
                )
                macro_ts[roi] = clean[:, 0]

            for roi in list(macro_ts_contra.keys()):
                clean = impute_nan_columns(
                    macro_ts_contra[roi][:, np.newaxis],
                    label=f"{subject}{run_tag} {label} contra {roi}"
                )
                macro_ts_contra[roi] = clean[:, 0]

            # ------------------------------------------------------
            # Step 2b — Collinearity diagnostic (NEW)
            #
            # Stack the full region pool used across all conditioning sets
            # for this subject/run/hemi (all loaded ipsi + contra regions)
            # and report the condition number of its covariance matrix.
            # This is a single upper-bound-ish proxy for how ill-conditioned
            # any individual conditioning set drawn from this pool could be
            # — it does not recompute a separate condition number for every
            # seed/target pair, which would be far more output for little
            # extra insight at this exploratory stage.
            # ------------------------------------------------------

            full_pool = (
                [macro_ts[r]        for r in loaded_rois] +
                [macro_ts_contra[r] for r in loaded_rois_contra]
            )
            X_pool = np.column_stack(full_pool)
            _condition_number_report(
                X_pool, subject, run_tag, label, n_time=X_pool.shape[0]
            )

            # ------------------------------------------------------
            # Step 3 — Partial correlations
            #
            # For every (seed, target) macro-region pair:
            #
            #   X = [seed | target | conditioning_macro-regions]
            #   C = partial_corr(X)
            #   result = C[0, 1]   ← seed ↔ target, all others partialled out
            #
            # Conditioning set = ALL OTHER macro-region timeseries, both
            # ipsilateral and contralateral, excluding seed and target.
            # This mirrors the task-free script where the bilateral MMP parcel
            # set is used as conditioning — inter-hemispheric confounds within
            # the task network are removed in the same way.
            #
            # Ipsilateral pairs:
            #   conditioning = remaining ipsi (excl. seed + target) + all contra
            #
            # Contralateral pairs:
            #   seed is ipsi, target is contra
            #   conditioning = remaining ipsi (excl. seed) + remaining contra (excl. target)
            #
            # Diagonal (seed == target) is left as NaN.
            # standardize='zscore_sample': safe on already-z-scored data
            # (near-no-op) and correct on non-z-scored individual runs.
            # ------------------------------------------------------

            n_rois            = len(clusters)
            partial_r_ipsi    = np.full((n_rois, n_rois), np.nan)
            partial_fz_ipsi   = np.full_like(partial_r_ipsi, np.nan)
            partial_r_contra  = np.full((n_rois, n_rois), np.nan)
            partial_fz_contra = np.full_like(partial_r_contra, np.nan)

            # cov_estimator is set by the ESTIMATOR_TAG toggle above (raw =
            # EmpiricalCovariance/unregularized, matching Nilearn's default;
            # ledoit-wolf = shrinkage-regularized). A fresh estimator instance
            # is created per fit_transform call inside Nilearn internally, so
            # reusing COV_ESTIMATOR across pairs here is safe.
            conn = ConnectivityMeasure(
                kind="partial correlation",
                cov_estimator=COV_ESTIMATOR,
                standardize=False,
            )

            for seed_name in loaded_rois:
                i_seed = clusters.index(seed_name)

                # --- Ipsilateral targets ---
                for tgt_name in loaded_rois:

                    if seed_name == tgt_name:
                        continue   # diagonal: leave as NaN

                    i_tgt = clusters.index(tgt_name)

                    # Conditioning: remaining ipsi (excl. seed + target) + all contra
                    ipsi_cond  = [macro_ts[r]        for r in loaded_rois        if r not in (seed_name, tgt_name)]
                    contra_cond = [macro_ts_contra[r] for r in loaded_rois_contra]
                    conditioning = ipsi_cond + contra_cond

                    if not conditioning:
                        print(f"  [{label}] ⚠️  No conditioning regions for {seed_name} → {tgt_name}; computing full correlation")

                    X = np.column_stack([macro_ts[seed_name], macro_ts[tgt_name]] + conditioning)

                    C = conn.fit_transform([X])[0]
                    r = C[0, 1]
                    partial_r_ipsi[i_seed,  i_tgt] = r
                    partial_fz_ipsi[i_seed, i_tgt] = np.arctanh(r)

                # --- Contralateral targets ---
                for tgt_name in loaded_rois_contra:

                    i_tgt = clusters.index(tgt_name)

                    # Conditioning: remaining ipsi (excl. seed) + remaining contra (excl. target)
                    ipsi_cond   = [macro_ts[r]        for r in loaded_rois        if r != seed_name]
                    contra_cond = [macro_ts_contra[r] for r in loaded_rois_contra if r != tgt_name]
                    conditioning = ipsi_cond + contra_cond

                    if not conditioning:
                        print(f"  [{label}] ⚠️  No conditioning regions for {seed_name} → contra {tgt_name}; computing full correlation")

                    X = np.column_stack([macro_ts[seed_name], macro_ts_contra[tgt_name]] + conditioning)

                    C = conn.fit_transform([X])[0]
                    r = C[0, 1]
                    partial_r_contra[i_seed,  i_tgt] = r
                    partial_fz_contra[i_seed, i_tgt] = np.arctanh(r)

            print(
                f"  [{label}] Partial corr: "
                f"{len(loaded_rois)} × {len(loaded_rois)} ipsi  |  "
                f"{len(loaded_rois)} × {len(loaded_rois_contra)} contra"
            )

            # ------------------------------------------------------
            # Step 4 — Map to full output grids
            #
            # Primary:   (n_clusters × n_clusters), ipsilateral only
            # Bilateral: (n_clusters × 2*n_clusters), [ipsi | contra]
            #            ipsi half always recoverable as [:, :n_clusters]
            #            column labels: {cluster}_ipsi / {cluster}_contra
            # ------------------------------------------------------

            n_cl = len(clusters)

            filled_r  = partial_r_ipsi.copy()    # already (n_cl × n_cl)
            filled_fz = partial_fz_ipsi.copy()

            filled_r_bilateral  = np.full((n_cl, 2 * n_cl), np.nan)
            filled_fz_bilateral = np.full_like(filled_r_bilateral, np.nan)

            # Ipsi half
            filled_r_bilateral[:,  :n_cl] = partial_r_ipsi
            filled_fz_bilateral[:, :n_cl] = partial_fz_ipsi

            # Contra half
            filled_r_bilateral[:,  n_cl:] = partial_r_contra
            filled_fz_bilateral[:, n_cl:] = partial_fz_contra

            clusters_bilateral = (
                [f"{c}_ipsi"  for c in clusters] +
                [f"{c}_contra" for c in clusters]
            )

            # ------------------------------------------------------
            # Step 5 — Save subject-level outputs
            # ------------------------------------------------------

            sub_out = (
                f"{main_data}/{subject}/91k/rest/corr/partial_corr"
                f"/by_hemi/task-constrained"
            )
            os.makedirs(sub_out, exist_ok=True)

            tag = label.lower()  # lh / rh

            # Build BIDS-like filename stem — ESTIMATOR_TAG appended at the
            # end so raw and ledoit-wolf outputs never overwrite each other.
            if run:
                stem = (
                    f"{subject}_task-rest_{run}"
                    f"_space-fsLR_den-91k_desc-fisher-z_{tag}_task-constrained"
                    f"_{ESTIMATOR_TAG}"
                )
            else:
                stem = (
                    f"{subject}"
                    f"_task-rest_space-fsLR_den-91k_desc-fisher-z_{tag}_task-constrained"
                    f"_{ESTIMATOR_TAG}"
                )

            # ---------- primary outputs ----------

            np.save(
                os.path.join(sub_out, f"{stem}.npy"),
                filled_fz,
            )

            pd.DataFrame(
                filled_fz,
                index=clusters,
                columns=clusters,
            ).to_csv(
                os.path.join(sub_out, f"{stem}.tsv"),
                sep="\t",
            )

            # ---------- bilateral outputs ----------

            np.save(
                os.path.join(sub_out, f"{stem}_bilateral.npy"),
                filled_fz_bilateral,
            )

            pd.DataFrame(
                filled_fz_bilateral,
                index=clusters,
                columns=clusters_bilateral,
            ).to_csv(
                os.path.join(sub_out, f"{stem}_bilateral.tsv"),
                sep="\t",
            )

            print(f"  [{label}] Saved to {sub_out}")

# ============================================================
# Collinearity diagnostic — final summary across all subjects
# ============================================================

print("\n" + "=" * 80)
print("COLLINEARITY DIAGNOSTIC SUMMARY — all subjects × runs × hemispheres")
print("=" * 80)

diag_df = pd.DataFrame(diagnostic_records)

if diag_df.empty:
    print("No diagnostic records collected (no hemispheres processed).")
else:
    # Sort worst-first (by RAW condition number) so problem cases are
    # immediately visible; Ledoit-Wolf column sits right next to it for
    # direct before/after comparison.
    diag_df_sorted = diag_df.sort_values("condition_number_raw", ascending=False)
    print(diag_df_sorted.to_string(index=False))

    print("\nSummary statistics — RAW (empirical) covariance:")
    print(diag_df["condition_number_raw"].describe().to_string())

    print("\nSummary statistics — Ledoit-Wolf shrunk covariance:")
    print(diag_df["condition_number_ledoit_wolf"].describe().to_string())

    n_total = len(diag_df)
    for col, flag_col, tag in (
        ("condition_number_raw",         "flag_raw",         "RAW"),
        ("condition_number_ledoit_wolf", "flag_ledoit_wolf", "Ledoit-Wolf"),
    ):
        n_high  = int((diag_df[flag_col] == "HIGH").sum())
        n_watch = int((diag_df[flag_col] == "WATCH").sum())
        n_ok    = int((diag_df[flag_col] == "OK").sum())
        print(
            f"\n{tag} flags: {n_high}/{n_total} HIGH (cond ≥ 100), "
            f"{n_watch}/{n_total} WATCH (30 ≤ cond < 100), "
            f"{n_ok}/{n_total} OK (cond < 30)"
        )

    median_raw = diag_df["condition_number_raw"].median()
    median_lw  = diag_df["condition_number_ledoit_wolf"].median()
    print(
        f"\nMedian condition number improved from {median_raw:.1f} (raw) "
        f"to {median_lw:.2f} (Ledoit-Wolf) "
        f"— {median_raw / median_lw:.0f}x reduction."
    )

    diag_out_dir = os.path.join(
        main_data, "group/91k/rest/partial_corr/diagnostics"
    )
    os.makedirs(diag_out_dir, exist_ok=True)
    # Note: this diagnostic computes BOTH raw and Ledoit-Wolf condition
    # numbers regardless of ESTIMATOR_TAG (see _condition_number_report),
    # so the file content is identical across runs — no need to tag it.
    diag_out_path = os.path.join(
        diag_out_dir, "collinearity_condition_numbers_task-constrained.csv"
    )
    diag_df.to_csv(diag_out_path, index=False)
    print(f"\nSaved diagnostic table: {diag_out_path}")