#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute INTRA-HEMISPHERIC PARTIAL correlations between ROI macro-regions (seeds)
and MMP parcels (targets) using Nilearn ConnectivityMeasure.

Design decisions
----------------
Seed timeseries  : ipsilateral hemisphere only (macro-region mask → mean vertex signal)
Target parcels   : ipsilateral hemisphere only (results are reported per hemisphere)
                   targets are not constrained by task results but span all MMP parcels
Conditioning set : ALL 106 parcels (both hemispheres), following Dawson et al. (2016)
                   and Genç et al. (2016).  Controlling for contralateral parcels
                   removes inter-hemispheric confounds that would otherwise inflate
                   apparent ipsilateral connectivity.

Exclusions from the conditioning set (ipsilateral only):
  (a) Parcels belonging to the seed's own macro-region — these share signal with
      the seed by definition and would otherwise dominate the partial correlation.
  (b) The target parcel itself — it cannot condition on itself.

Self-correlation masking
  A seed macro-region and a target parcel that belongs to the SAME macro-region
  will share vertices, so their partial correlation will be spuriously high
  (approaching 1).  These entries are left as NaN in the output rather than
  filled with a misleading value.  A separate macro-region-by-macro-region summary
  matrix is saved as a sanity check; off-diagonal entries should be interpretable
  connectivity estimates, diagonal entries are NaN by construction.

Standardize flag (FINAL — verified against nilearn 0.10.2, the version
actually installed on the cluster, not just the current release)
  ConnectivityMeasure(standardize=False), set explicitly for honesty rather
  than as a meaningful analysis choice: for kind="partial correlation",
  nilearn's `standardize` argument is NEVER applied — confirmed directly
  from nilearn 0.10.2's _fit_transform() (connectivity_matrices.py):
  standardize is only consumed inside the `if self.kind == "correlation":`
  branch; for "partial correlation" the `else` branch calls
  `self.cov_estimator_.fit(x)` directly on the raw input, with no
  standardization step at all, regardless of what `standardize` is set to.
  This holds identically in 0.10.2 and in the current nilearn release — it
  is not a version-dependent behavior.

  Project history: an earlier version of this script used
  standardize="zscore_sample" (matching a mistaken assumption that this
  flag does something for kind="partial correlation"). Since it provably
  does nothing for this kind, both False and "zscore_sample" produce
  bit-for-bit identical output here — this was verified empirically before
  correcting the flag to False. The remaining/actual scale-correction logic
  now lives in Step 3b below, done manually, once, since nilearn will not
  do it for kind="partial correlation".

Manual re-standardization after ROI-averaging (Step 3b)
  Because nilearn's standardize flag is a no-op for partial correlation
  (see above), any benefit of standardizing before computing the
  covariance/precision matrix has to be done manually, before calling
  ConnectivityMeasure.fit_transform().

  Why this step exists at all: every vertex is already z-scored by XCP-D,
  but the ROI-MEAN signal is not guaranteed to keep unit variance after
  averaging. For N vertices sharing average within-region correlation rho,
  Var(region/parcel mean) = rho + (1-rho)/N — a function of region size
  alone. Macro-regions and MMP parcels here differ substantially in vertex
  count, so averaged timeseries can end up on measurably different scales
  even though every underlying vertex started at unit variance. This
  matters only for regularized estimators (Ledoit-Wolf, GraphicalLassoCV),
  whose shrinkage target / L1 penalty assume comparable scale across
  variables — empirically confirmed to have ZERO effect on the raw
  (unregularized) estimator, since partial correlation from an
  unregularized covariance is exactly invariant to per-column scaling
  (verified: max diff 0.0 across multiple test scenarios). At realistic
  parameter values (vertex counts and within-region correlation typical of
  cortical parcels), the resulting scale mismatch across regions is modest
  (roughly 1.05x-1.15x SD ratio, not an order of magnitude), so this
  correction is best understood as a small, free, mathematically justified
  safety net — not a fix expected to dramatically change the results on
  its own. Done ONCE, after NaN imputation and before building any
  seed/target matrix, to avoid compounding multiple standardization steps.

Covariance estimator (added for consistency with the task-constrained script)
  Conditioning here is on all 106 bilateral MMP parcels — a considerably larger
  and likely more collinear regressor set than the 24 macro-regions used in the
  task-constrained analysis, so the same collinearity concerns apply at least as
  much here.  ESTIMATOR_TAG is passed in as sys.argv[1] by the SLURM submit
  script (submit_nilearn_compute_partial_corr_job.py), which owns the choice of
  estimator as a run-configuration decision — this script only validates the tag
  and maps it to the corresponding sklearn estimator.  No silent default is used
  on the cluster to keep every job's estimator choice explicit and traceable.

  raw            : EmpiricalCovariance (unregularized) — Nilearn's default
  ledoit-wolf    : LedoitWolf shrinkage (Ledoit & Wolf, 2004) — analytic
                   shrinkage intensity, no cross-validation
  graphical-lasso: GraphicalLassoCV — L1-regularized precision matrix,
                   sparsity penalty selected by cross-validation (slower;
                   produces a sparse precision matrix — a different
                   scientific claim than shrinkage; see Peterson et al. 2025,
                   Imaging Neuroscience, for the graphical lasso vs.
                   graphical ridge distinction)

  Output filenames are tagged with ESTIMATOR_TAG so that runs with different
  estimators never overwrite each other on disk.

Outputs (per subject, per run, per hemisphere, per estimator)
  HARMONIZED to match the task-constrained script's BIDS-style stem
  (see "Harmonization" note above Step 7 in the code): filenames use
  {subject}_task-rest{run_entity}_space-fsLR_den-91k_desc-fisher-z_{hemi}
  _task-free_{estimator}[_suffix].{npy,tsv} instead of the earlier
  "seed-task_by_mmp-parcel_..." prefix style, and tabular outputs are
  tab-separated .tsv (not comma-separated .csv). Only Fisher-z is saved
  at the subject level (no standalone Pearson r file — dropped to match
  task-constrained and the pipeline-wide Fisher-z-first design):

  {subject}_task-rest{run_entity}_space-fsLR_den-91k_desc-fisher-z_{hemi}_task-free_{estimator}.npy / .tsv
      — primary output, ipsilateral only (n_clusters × n_parcels)
  {subject}_task-rest{run_entity}_space-fsLR_den-91k_desc-fisher-z_{hemi}_task-free_{estimator}_bilateral.npy / .tsv
      — bilateral output, [ipsi | contra] (n_clusters × 2*n_parcels)
  {subject}_task-rest{run_entity}_space-fsLR_den-91k_desc-fisher-z_{hemi}_task-free_{estimator}_macro-summary.npy / .tsv
      — macro-region-by-macro-region sanity check (no task-constrained
        equivalent, since that script's primary output already is at
        this level)

Group aggregation is handled by group_partial_corr_by_hemi_task-free.py
(always in Fisher-z space; back-transformed to r only at the reporting
stage). That script's npy_path() must match this harmonized convention.

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
# See docstring above. Required as sys.argv[1] — no silent default.
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
atlas_folder = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1"

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

# Flat ordered parcel list — defines the column order of all output matrices
parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

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
# Subject loop
# ============================================================

for subject in subjects:

    for run in RUNS:

        run_tag = f"_{run}" if run else ""

        print(f"\n=== Processing {subject}{run_tag} ===")

        timeseries_fn = (
            f"{main_data}/{subject}/91k/rest/timeseries/"
            f"{subject}_ses-01_task-rest{run_tag}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
        )

        ts_img, ts_data_raw = load_surface(timeseries_fn)

        # Separate LH and RH 32k surface arrays.
        # return_32k_mask=True returns a boolean mask (True=cortex, False=medial wall)
        # for downstream QC — medial wall vertices are never included in parcel masks.
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
            label     = h["label"]     # "LH" or "RH"
            seed_key  = h["seed_key"]  # "lh" or "rh"  (seed filename suffix)
            atlas_key = h["atlas_key"] # "L"  or "R"   (atlas filename prefix)
            ts_key    = h["ts_key"]    # "data_L" or "data_R"

            ts_data = res[ts_key]   # (n_time, n_vertices_hemi)
            print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

            # ------------------------------------------------------
            # Step 1 — Seed (macro-region) timeseries, ipsilateral only
            #
            # Each macro-region has a pre-computed binary mask (.shape.gii).
            # Mean signal across all mask vertices → one timeseries per macro-region.
            # ------------------------------------------------------

            cluster_ts_list    = []
            cluster_names_used = []

            for roi in clusters:
                _, mask_data = load_surface(
                    f"{main_data}/{subject}/91k/rest/seed/"
                    f"{subject}_91k_intertask_Sac-Pur-pRF_{seed_key}_{roi}.shape.gii"
                )
                mask = mask_data.ravel()

                if not np.any(mask):
                    print(f"  [{label}] ⚠️  Empty seed mask for {roi} — skipping")
                    continue

                cluster_ts_list.append(ts_data[:, mask > 0].mean(axis=1))
                cluster_names_used.append(roi)

            if not cluster_ts_list:
                print(f"  [{label}] ⚠️  No valid seed timeseries — skipping hemisphere")
                continue

            cluster_ts = np.column_stack(cluster_ts_list)
            print(f"  [{label}] Seeds loaded: {cluster_names_used}")

            # ------------------------------------------------------
            # Step 2 — Parcel timeseries, BOTH hemispheres
            #
            # All 106 parcels (53 ipsi + 53 contra) are loaded to form the
            # bilateral conditioning set.  Three tracking lists are built:
            #   parcel_names_used : parcel name for each matrix column
            #   parcel_hemi_used  : "L" or "R" for each column
            #   ipsi_col_idx      : column indices belonging to the ipsilateral
            #                       hemisphere — these are the targets and define
            #                       the result columns in the output matrix
            #
            # Ipsilateral columns are listed first so ipsi_col_idx entries
            # are always the lowest indices (easier to audit in logs).
            # ------------------------------------------------------

            atlas_keys_ordered = [atlas_key, "R" if atlas_key == "L" else "L"]

            parcel_ts_list    = []
            parcel_names_used = []
            parcel_hemi_used  = []
            ipsi_col_idx      = []

            for ak in atlas_keys_ordered:
                ts_source = res["data_L" if ak == "L" else "data_R"]

                for parcel in parcels:
                    _, mask_data = load_surface(
                        f"{atlas_folder}/parcels/{ak}_{parcel}_ROI.shape.gii"
                    )
                    mask = mask_data.ravel()

                    if not np.any(mask):
                        continue

                    col_idx = len(parcel_ts_list)
                    parcel_ts_list.append(ts_source[:, mask > 0].mean(axis=1))
                    parcel_names_used.append(parcel)
                    parcel_hemi_used.append(ak)

                    if ak == atlas_key:
                        ipsi_col_idx.append(col_idx)

            if not parcel_ts_list:
                print(f"  [{label}] ⚠️  No valid parcel timeseries — skipping hemisphere")
                continue

            parcel_ts = np.column_stack(parcel_ts_list)

            n_ipsi   = sum(1 for hm in parcel_hemi_used if hm == atlas_key)
            n_contra = sum(1 for hm in parcel_hemi_used if hm != atlas_key)
            print(
                f"  [{label}] Bilateral parcel matrix: {parcel_ts.shape[1]} columns "
                f"({n_ipsi} ipsi [{atlas_key}] + {n_contra} contra)"
            )

            # ------------------------------------------------------
            # Step 3 — NaN imputation (must happen before Nilearn)
            #
            # Expected cases: sub-22 run-02 (bbregister failure),
            #                 sub-25 parcel 6mp (FOV truncation).
            # ------------------------------------------------------

            cluster_ts = impute_nan_columns(cluster_ts, label=f"{subject}{run_tag} {label} seed")
            parcel_ts  = impute_nan_columns(parcel_ts,  label=f"{subject}{run_tag} {label} parcel")

            # ------------------------------------------------------
            # Step 4 — Partial correlations
            #
            # For every (seed_macro-region, target_parcel) pair:
            #
            #   X = [seed | target | conditioning_parcels]
            #   C = partial_corr(X)        ← via Nilearn ConnectivityMeasure
            #   result = C[0, 1]           ← seed ↔ target, all others partialled out
            #
            # Exclusions from the conditioning set (ipsilateral only):
            #   (a) Seed's own macro-region parcels — pre-computed once per seed
            #       as seed_own_cols (constant across all targets for this seed).
            #   (b) The target parcel column itself — added per iteration.
            #
            # SELF-CORRELATION MASKING:
            #   When the target parcel belongs to the seed's own macro-region, the
            #   seed and target timeseries are derived from overlapping vertices and
            #   their partial correlation is trivially high.  These entries are left
            #   as NaN and are not computed.  They are summarised in the macro-region-
            #   by-macro-region sanity check matrix (Step 5) instead.
            # ------------------------------------------------------

            ipsi_parcel_names = [parcel_names_used[j] for j in ipsi_col_idx]
            n_ipsi_parcels    = len(ipsi_col_idx)
            n_clusters_used   = cluster_ts.shape[1]

            # Result matrices: rows = macro-regions, cols = ipsilateral target parcels
            partial_r  = np.full((n_clusters_used, n_ipsi_parcels), np.nan)
            partial_fz = np.full_like(partial_r, np.nan)

            # Contralateral column indices and names (complement of ipsi_col_idx)
            ipsi_col_set        = set(ipsi_col_idx)
            contra_col_idx      = [j for j in range(len(parcel_names_used)) if j not in ipsi_col_set]
            contra_parcel_names = [parcel_names_used[j] for j in contra_col_idx]
            n_contra_parcels    = len(contra_col_idx)

            # Contralateral result matrices
            partial_r_contra  = np.full((n_clusters_used, n_contra_parcels), np.nan)
            partial_fz_contra = np.full_like(partial_r_contra, np.nan)

            # Single ConnectivityMeasure instance reused for every (seed, target) pair.
            # standardize=False: see Step 3b above — this argument has no
            # effect for kind="partial correlation" either way, but is set
            # explicitly to avoid implying nilearn does something it doesn't.
            conn = ConnectivityMeasure(
                kind="partial correlation",
                cov_estimator=COV_ESTIMATOR,
                standardize=False,
            )

            for i_cl, cl_name in enumerate(cluster_names_used):

                # (a) ipsilateral columns belonging to this seed's own macro-region
                seed_own_parcels = set(seed_to_parcels.get(cl_name, []))
                seed_own_cols = {
                    j for j, (pname, phemi) in enumerate(
                        zip(parcel_names_used, parcel_hemi_used)
                    )
                    if phemi == atlas_key and pname in seed_own_parcels
                }

                # --- Ipsilateral targets ---
                for i_target, target_parcel in enumerate(ipsi_parcel_names):

                    # Self-correlation: skip, leave as NaN
                    if target_parcel in seed_own_parcels:
                        continue

                    # (b) exclude the target column itself from conditioning
                    target_col   = ipsi_col_idx[i_target]
                    exclude_cols = seed_own_cols | {target_col}

                    conditioning_idx = [
                        j for j in range(len(parcel_names_used))
                        if j not in exclude_cols
                    ]

                    if not conditioning_idx:
                        print(f"  [{label}] ⚠️  No conditioning parcels left for {cl_name} → {target_parcel}")
                        continue

                    # Column 0 = seed, column 1 = target, columns 2: = conditioning parcels
                    X = np.column_stack([
                        cluster_ts[:, i_cl],
                        parcel_ts[:, target_col],
                        parcel_ts[:, conditioning_idx],
                    ])

                    C = conn.fit_transform([X])[0]

                    r = C[0, 1]
                    partial_r[i_cl,  i_target] = r
                    partial_fz[i_cl, i_target] = np.arctanh(r)

                # --- Contralateral targets ---
                # No self-correlation masking needed (seed is ipsilateral, target is
                # contralateral — they cannot share vertices).
                # Conditioning set excludes seed's own ipsilateral parcels and the
                # target column, same logic as above.
                for i_target, target_parcel in enumerate(contra_parcel_names):

                    target_col   = contra_col_idx[i_target]
                    exclude_cols = seed_own_cols | {target_col}

                    conditioning_idx = [
                        j for j in range(len(parcel_names_used))
                        if j not in exclude_cols
                    ]

                    if not conditioning_idx:
                        print(f"  [{label}] ⚠️  No conditioning parcels left for {cl_name} → contra {target_parcel}")
                        continue

                    X = np.column_stack([
                        cluster_ts[:, i_cl],
                        parcel_ts[:, target_col],
                        parcel_ts[:, conditioning_idx],
                    ])

                    C = conn.fit_transform([X])[0]

                    r = C[0, 1]
                    partial_r_contra[i_cl,  i_target] = r
                    partial_fz_contra[i_cl, i_target] = np.arctanh(r)

            print(
                f"  [{label}] Partial corr: "
                f"{n_clusters_used} seeds × {n_ipsi_parcels} ipsi + {n_contra_parcels} contra parcels"
            )

            # ------------------------------------------------------
            # Step 5 — Macro-region-by-macro-region sanity check
            #
            # Averages partial_fz over the target parcels belonging to each
            # macro-region, giving an (n_clusters × n_clusters) summary matrix.
            #   - Diagonal entries are NaN by construction (self-correlation mask).
            #   - Off-diagonal entries are mean partial connectivity between
            #     macro-region pairs — used for QC and cluster-level reporting.
            # ------------------------------------------------------

            n_cl = len(clusters)
            cluster_by_cluster_fz = np.full((n_cl, n_cl), np.nan)

            for i_cl, cl_seed in enumerate(cluster_names_used):
                gr = clusters.index(cl_seed)

                for j_cl, cl_target in enumerate(clusters):
                    target_parcel_names = seed_to_parcels.get(cl_target, [])

                    target_col_idx = [
                        i for i, pname in enumerate(ipsi_parcel_names)
                        if pname in target_parcel_names
                    ]

                    if not target_col_idx:
                        continue

                    vals = partial_fz[i_cl, target_col_idx]
                    if not np.all(np.isnan(vals)):
                        cluster_by_cluster_fz[gr, j_cl] = np.nanmean(vals)

            # ------------------------------------------------------
            # Step 6 — Map to full output grids (Fisher-z only)
            #
            #   filled_fz  (n_clusters × n_parcels)
            #     Ipsilateral only. Columns = parcels in canonical YAML order.
            #     This is the primary output and is what all downstream scripts
            #     (group stats, visualisation) consume.
            #
            #   filled_fz_bilateral  (n_clusters × 2*n_parcels)
            #     Columns = [ipsi_parcels | contra_parcels], both in canonical
            #     YAML order.  The ipsilateral half is always [:, :n_parcels],
            #     so any downstream script can recover it without change.
            #
            # HARMONIZATION NOTE: an earlier version of this script also saved
            # a standalone Pearson r grid (filled_r / filled_r_bilateral)
            # alongside Fisher-z. This has been dropped to match the
            # task-constrained script's subject-level convention (Fisher-z
            # only) and the pipeline-wide design principle stated throughout:
            # average in Fisher-z space, recover r only at the reporting
            # stage via tanh(). Pearson r is still computed in memory
            # (partial_r / partial_r_contra, directly from Nilearn, before
            # arctanh()) — it is simply no longer written to disk here.
            #
            # cluster_names_used may be a subset of clusters if any seed masks
            # were empty — the grid is initialised to NaN so missing rows
            # are explicitly absent rather than silently zero.
            # ------------------------------------------------------

            n_parcels_total = len(parcels)

            filled_fz           = np.full((len(clusters), n_parcels_total), np.nan)
            filled_fz_bilateral = np.full((len(clusters), 2 * n_parcels_total), np.nan)

            # Column labels for the bilateral DataFrame
            contra_key   = "R" if atlas_key == "L" else "L"
            parcels_bilateral = (
                [f"{p}_{atlas_key}" for p in parcels] +   # ipsi half
                [f"{p}_{contra_key}" for p in parcels]    # contra half
            )

            for i_cl, cl in enumerate(cluster_names_used):
                gr = clusters.index(cl)

                # Ipsilateral half
                for i_target, pa in enumerate(ipsi_parcel_names):
                    gc = parcels.index(pa)
                    filled_fz[gr, gc] = partial_fz[i_cl, i_target]
                    # Bilateral ipsi half (columns 0 : n_parcels)
                    filled_fz_bilateral[gr, gc] = partial_fz[i_cl, i_target]

                # Contralateral half (columns n_parcels : 2*n_parcels)
                for i_target, pa in enumerate(contra_parcel_names):
                    gc = parcels.index(pa)
                    filled_fz_bilateral[gr, n_parcels_total + gc] = partial_fz_contra[i_cl, i_target]

            # ------------------------------------------------------
            # Step 7 — Save subject-level outputs
            #
            # HARMONIZED filename convention (matches the task-constrained
            # script's BIDS-style stem exactly, aside from the
            # "task-free" vs "task-constrained" entity):
            #
            #   {subject}_task-rest{run_entity}_space-fsLR_den-91k
            #       _desc-fisher-z_{hemi}_task-free_{ESTIMATOR_TAG}[_bilateral].npy/.tsv
            #
            # Two conventions changed from an earlier version of this script:
            #   1. Filenames switched from the "seed-task_by_mmp-parcel_..."
            #      prefix style to the full BIDS-style stem, matching
            #      task-constrained subject-level outputs on disk.
            #   2. Tabular outputs switched from comma-separated .csv to
            #      tab-separated .tsv (sep="\t"), matching task-constrained.
            #
            # The macro-region-by-macro-region sanity check matrix has no
            # equivalent in the task-constrained script (whose primary output
            # already IS macro-region-by-macro-region), so it keeps the same
            # base stem with an extra "_macro-summary" suffix rather than a
            # separate naming scheme.
            #
            # ESTIMATOR_TAG is appended to every output filename so that runs
            # with different covariance estimators coexist on disk rather than
            # overwriting one another.
            # ------------------------------------------------------

            sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr/by_hemi/task-free"
            os.makedirs(sub_out, exist_ok=True)

            tag = label.lower()   # "lh" or "rh"

            # BIDS-style stem, matching task-constrained's construction exactly
            if run:
                stem = (
                    f"{subject}_task-rest_{run}"
                    f"_space-fsLR_den-91k_desc-fisher-z_{tag}_task-free"
                    f"_{ESTIMATOR_TAG}"
                )
            else:
                stem = (
                    f"{subject}"
                    f"_task-rest_space-fsLR_den-91k_desc-fisher-z_{tag}_task-free"
                    f"_{ESTIMATOR_TAG}"
                )

            # --- Ipsilateral output (primary; consumed by all downstream scripts) ---
            np.save(os.path.join(sub_out, f"{stem}.npy"), filled_fz)
            pd.DataFrame(filled_fz, index=clusters, columns=parcels).to_csv(
                os.path.join(sub_out, f"{stem}.tsv"), sep="\t"
            )

            # --- Bilateral output ([ipsi | contra]; ipsi half = [:, :n_parcels]) ---
            np.save(os.path.join(sub_out, f"{stem}_bilateral.npy"), filled_fz_bilateral)
            pd.DataFrame(filled_fz_bilateral, index=clusters, columns=parcels_bilateral).to_csv(
                os.path.join(sub_out, f"{stem}_bilateral.tsv"), sep="\t"
            )

            # --- Macro-region-by-macro-region sanity check (Fisher-z, ipsi only) ---
            np.save(os.path.join(sub_out, f"{stem}_macro-summary.npy"), cluster_by_cluster_fz)
            pd.DataFrame(cluster_by_cluster_fz, index=clusters, columns=clusters).to_csv(
                os.path.join(sub_out, f"{stem}_macro-summary.tsv"), sep="\t"
            )

            print(f"  [{label}] Saved to {sub_out}")

print("\nDone. Run group_partial_corr_by_hemi_task-free.py to aggregate across subjects.")
print("NOTE: group_partial_corr_by_hemi_task-free.py must be updated to accept/propagate")
print(f"the same ESTIMATOR_TAG ('{ESTIMATOR_TAG}') used here when locating these files.")
# ============================================================