"""
-----------------------------------------------------------------------------------------
make_corr_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute per-vertex pRF derivative correlations (position, size, R2) between two SacLoc /
pRF conditions, ROI-by-ROI, for either a single subject or a group of subjects.
The two conditions can be:
    - pmf vs pmf   (task-SacLoc, e.g. pmf-css-rdm  vs  pmf-css-odm)
    - pmf vs prf   (e.g. pmf-css-rdm  vs  prf-css,  or  pmf-css-odm  vs  prf-css)

For an individual subject (sys.argv[3] does not contain "group"):
- Per-vertex derivative maps are loaded from the surface (.func.gii or .dtseries.nii) for
  each condition, split into ROIs (rois-group-mmp), and merged on shared vertices (inner
  join on roi/hemi/num_vert).
- The merged, un-thresholded per-vertex table is saved as a TSV (useful for QC / re-use).
- For each parameter and ROI, vertices below RSQ_THRESHOLD in either condition are
  dropped, then a weighted 2D KDE is computed on a FIXED grid_size x grid_size grid
  (grid range = corr_param_ranges below -- fixed on purpose, see note below), weighted by
  the mean R2 of both conditions.
- A weighted Deming (orthogonal) regression is computed per ROI, weighted the same way,
  along with Pearson r2, p-value, n_vertex and mean R2 weight.

For a group (sys.argv[3] contains "group", e.g. "group-all"):
- This does NOT reload any surface data. It averages the per-subject TSVs that this same
  script must already have produced by being run once per individual subject first.
- KDE density grids are averaged cell-by-cell across subjects (every subject counts
  equally, regardless of how many vertices it contributed) -- this only works because the
  grid is identical for every subject (see note below).
- Regression stats (slope, intercept, r2_pearson) are averaged across subjects, with
  2.5 / 97.5 percentile bounds across subjects reported as a between-subject CI.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: subject (e.g. sub-11) OR group name (e.g. group-all, defined in settings.yml)
sys.argv[4]: condition A (e.g. pmf-css-rdm)
sys.argv[5]: condition B (e.g. pmf-css-odm, or prf-css for a pmf-vs-prf comparison)
-----------------------------------------------------------------------------------------
Output(s), all under {main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format}/pmf/tsv :
- {subject}_{cond_a}_vs_{cond_b}_deriv.tsv
      merged per-vertex table for the two conditions (individual subjects only)
- {subject}_{cond_a}_vs_{cond_b}_{corr_param}-corr.tsv
      columns: roi, x_grid, y_grid, density
- {subject}_{cond_a}_vs_{cond_b}_{corr_param}-regression.tsv
      columns: roi, slope, intercept, r2_pearson, p_value, n_vertex, mean_r2_weight
      (+ *_ci_lower / *_ci_upper for group runs)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python make_corr_tsv.py /scratch/mszinte/data RetinoMaps sub-11 pmf-css-rdm pmf-css-odm
python make_corr_tsv.py /scratch/mszinte/data RetinoMaps group-all pmf-css-rdm pmf-css-odm
-----------------------------------------------------------------------------------------
Adapted from Martin Szinte's make_corr_tsv.py / make_corr_fig.py (amblyo7T_prf) 
-----------------------------------------------------------------------------------------
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, pearsonr

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import load_surface
from pycortex_utils import get_rois, set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
cond_a = sys.argv[4]  # e.g. "pmf-css-rdm"
cond_b = sys.argv[5]  # e.g. "pmf-css-odm" or "prf-css"

RSQ_THRESHOLD = 0.1
GRID_SIZE = 100
rois_method_format = "rois-group-mmp"

TASK_NAME_BY_ANALYSIS = {"pmf": "SacLoc", "prf": "pRF"}

rois = ["V1", "V2", "V3", "V3AB", "LO", "VO", "hMT+", "iIPS", "sIPS", "iPCS", "sPCS", "mPCS"]

corr_params = ["prf_x", "prf_y", "prf_ecc", "prf_size", "prf_rsq"]
corr_param_ranges = {
    "prf_x": [-20, 20],
    "prf_y": [-20, 20],
    "prf_ecc": [0, 20],
    "prf_size": [0, 10],
    "prf_rsq": [0, 0.8],
}

# "pmf-css-rdm" -> "PMF" ; "prf-css" -> "prf" / "PRF" (only one prf condition)
analysis_type_a, analysis_type_b = cond_a.split("-")[0], cond_b.split("-")[0]
label_a = cond_a.split("-")[-1].upper() if analysis_type_a == "pmf" else "PRF"
label_b = cond_b.split("-")[-1].upper() if analysis_type_b == "pmf" else "PRF"

# Load settings -- pmf-vs-prf comparisons need both analysis-specific yml files
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
analysis_info_by_type = {}
for atype in sorted({analysis_type_a, analysis_type_b}):
    analysis_settings_path = os.path.join(base_dir, project_dir, f"{atype}-analysis.yml")
    settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
    analysis_info_by_type[atype] = settings[0]

analysis_info = analysis_info_by_type[analysis_type_a]
formats = analysis_info["formats"]
maps_names_css = analysis_info["maps_names_css"]

# ── weighted Deming (orthogonal) regression, symmetric when x and y are swapped ─────────
def weighted_deming_regression(x, y, weights=None):
    x, y = np.array(x), np.array(y)
    weights = np.ones(len(x)) if weights is None else np.array(weights)
    weights = weights / weights.sum()
    x_mean, y_mean = np.sum(weights * x), np.sum(weights * y)
    sxx = np.sum(weights * (x - x_mean) ** 2)
    syy = np.sum(weights * (y - y_mean) ** 2)
    sxy = np.sum(weights * (x - x_mean) * (y - y_mean))
    slope = (syy - sxx + np.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)
    intercept = y_mean - slope * x_mean
    return slope, intercept


# ═════════════════════════════════════════════════════════════════════════════════════
# INDIVIDUAL SUBJECT
# ═════════════════════════════════════════════════════════════════════════════════════
if "group" not in subject:

    for format_ in formats:
        print(f"\n=== format: {format_} ===")
        tsv_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/tsv"
        os.makedirs(tsv_dir, exist_ok=True)
        pycortex_cortex_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/cortex"
        set_pycortex_config_file(pycortex_cortex_dir)
        # 170k (CIFTI) surfaces are registered to a template pycortex subject, not the
        # individual subject's own fsnative surface
        pycortex_subject = analysis_info["pycortex_subject_template"] if format_ == "170k" else subject

        # ── 1. load per-vertex data for the two conditions, split by ROI ────────────
        dfs = {}
        for cond, atype in ((cond_a, analysis_type_a), (cond_b, analysis_type_b)):
            info = analysis_info_by_type[atype]
            task_name = TASK_NAME_BY_ANALYSIS[atype]
            preproc_prep, filtering, normalization = info["preproc_prep"], info["filtering"], info["normalization"]
            avg_method = info["avg_methods"][0]
            maps_names_css = info["maps_names_css_loo"] if "loo" in avg_method else info["maps_names_css"]
            if format_ == "170k":
                # 170k (CIFTI) derivatives are one combined file, no hemi split
                deriv_input = (
                    f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{atype}/"
                    f"prf_derivatives/{subject}_task-{task_name}_{preproc_prep}_{filtering}_"
                    f"{normalization}_{avg_method}_{cond}_deriv.dtseries.nii"
                )
            else:
                deriv_input = {
                    hemi: (
                        f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{atype}/"
                        f"prf_derivatives/{subject}_task-{task_name}_{hemi}_{preproc_prep}_{filtering}_"
                        f"{normalization}_{avg_method}_{cond}_deriv.func.gii"
                    )
                    for hemi in ("hemi-L", "hemi-R")
                }
            print(f"loading {cond} for {subject} (avg_method={avg_method})")

            # turn the derivative surface map(s) into a per-ROI, per-vertex dataframe
            df_cond = pd.DataFrame()
            items = deriv_input.items() if isinstance(deriv_input, dict) else [(None, deriv_input)]
            for hemi, deriv_fn in items:
                _, deriv_mat = load_surface(deriv_fn)  # (n_params, n_vertices)
                assert deriv_mat.shape[0] == len(maps_names_css), \
                    f"{deriv_fn} has {deriv_mat.shape[0]} frames, expected {len(maps_names_css)} — check loo_r2"
                get_rois_kwargs = dict(subject=pycortex_subject, surf_format=format_,
                                        rois_type=rois_method_format, mask=True, rois=rois)
                if hemi is not None:
                    get_rois_kwargs["hemis"] = hemi
                roi_verts = get_rois(**get_rois_kwargs)
                for roi in roi_verts.keys():
                    data_dict = {col: deriv_mat[i, roi_verts[roi]] for i, col in enumerate(maps_names_css)}
                    data_dict["roi"] = roi
                    data_dict["hemi"] = hemi if hemi is not None else "both"
                    data_dict["num_vert"] = np.where(roi_verts[roi])[0]
                    df_cond = pd.concat([df_cond, pd.DataFrame(data_dict)], ignore_index=True)
            dfs[cond] = df_cond

        # ── 2. merge the two conditions on shared vertices ───────────────────────────
        df_a, df_b = dfs[cond_a], dfs[cond_b]
        data = pd.merge(df_a, df_b, on=["roi", "hemi", "num_vert"], suffixes=(f"_{label_a}", f"_{label_b}"))
        print(f"Merged {len(data)} vertices across {data['roi'].nunique()} rois ({label_a} vs {label_b})")

        # diagnostic only -- does NOT change the (fixed) grid used below
        for p in ("prf_x", "prf_y", "prf_ecc"):
            bound = corr_param_ranges[p][1]
            for lbl in (label_a, label_b):
                v = data[f"{p}_{lbl}"]
                n_railed = (v.abs() >= bound - 1e-6).sum()
                if n_railed > 1:
                    print(f"  NOTE: {p}_{lbl} has {n_railed} vertices at |value|~={bound:.3f} "
                          f"— check whether this is a grid-search bound (railed fits)")

        deriv_fn = f"{tsv_dir}/{subject}_{cond_a}_vs_{cond_b}_deriv.tsv"
        data.to_csv(deriv_fn, sep="\t", na_rep="NaN", index=False)
        print(f"Saving: {deriv_fn}")

        # ── 3. per-parameter, per-ROI weighted KDE + weighted Deming regression ──────
        for corr_param in corr_params:
            param_range = corr_param_ranges[corr_param]
            x_grid = np.linspace(param_range[0], param_range[1], GRID_SIZE)
            y_grid = np.linspace(param_range[0], param_range[1], GRID_SIZE)
            xx, yy = np.meshgrid(x_grid, y_grid)
            grid_points = np.vstack([xx.ravel(), yy.ravel()])

            kde_rows, reg_rows = [], []
            for roi in rois:
                df_roi = data.loc[data["roi"] == roi]
                x_vals = df_roi[f"{corr_param}_{label_a}"].values
                y_vals = df_roi[f"{corr_param}_{label_b}"].values
                r2_a = df_roi[f"prf_rsq_{label_a}"].values
                r2_b = df_roi[f"prf_rsq_{label_b}"].values
                r2_combined = (r2_a + r2_b) / 2
                # exclude NaNs AND sub-threshold/negative-rsq vertices — gaussian_kde's
                # weights must be >= 0, and low-rsq vertices are noise anyway
                mask = (
                    ~(np.isnan(x_vals) | np.isnan(y_vals) | np.isnan(r2_combined))
                    & (r2_a > RSQ_THRESHOLD) & (r2_b > RSQ_THRESHOLD)
                )
                x_vals, y_vals, r2_combined = x_vals[mask], y_vals[mask], r2_combined[mask]
                n_vertex = len(x_vals)
                mean_r2_weight = np.nanmean(r2_combined) if n_vertex > 0 else np.nan

                if n_vertex > 1:
                    r2_norm = r2_combined / r2_combined.sum()
                    try:
                        kde = gaussian_kde(np.vstack([x_vals, y_vals]), weights=r2_norm)
                        density = kde(grid_points).reshape(GRID_SIZE, GRID_SIZE)
                    except np.linalg.LinAlgError:
                        density = np.full((GRID_SIZE, GRID_SIZE), np.nan)
                else:
                    r2_norm, density = np.array([]), np.full((GRID_SIZE, GRID_SIZE), np.nan)

                kde_rows.append(pd.DataFrame({
                    "roi": roi, "x_grid": xx.ravel(), "y_grid": yy.ravel(), "density": density.ravel(),
                }))

                if n_vertex > 2:
                    slope, intercept = weighted_deming_regression(x_vals, y_vals, r2_norm)
                    r_val, p_val = pearsonr(x_vals, y_vals)
                    r2_pearson = r_val ** 2
                else:
                    slope = intercept = r2_pearson = p_val = np.nan
                reg_rows.append({
                    "roi": roi, "slope": slope, "intercept": intercept, "r2_pearson": r2_pearson,
                    "p_value": p_val, "n_vertex": n_vertex, "mean_r2_weight": mean_r2_weight,
                })

            df_kde = pd.concat(kde_rows, ignore_index=True)
            df_reg = pd.DataFrame(reg_rows)

            kde_fn = f"{tsv_dir}/{subject}_{cond_a}_vs_{cond_b}_{corr_param}-corr.tsv"
            df_kde.to_csv(kde_fn, sep="\t", na_rep="NaN", index=False)
            print(f"Saving tsv: {kde_fn}")

            reg_fn = f"{tsv_dir}/{subject}_{cond_a}_vs_{cond_b}_{corr_param}-regression.tsv"
            df_reg.to_csv(reg_fn, sep="\t", na_rep="NaN", index=False)
            print(f"Saving tsv: {reg_fn}")

# ═════════════════════════════════════════════════════════════════════════════════════
# GROUP -- averages per-subject TSVs already produced by the block above
# ═════════════════════════════════════════════════════════════════════════════════════
else:
    subjects_to_group = analysis_info["subjects"]

    for format_ in formats:
        print(f"\n=== format: {format_} — group: {subject} ({len(subjects_to_group)} subjects) ===")
        tsv_dir_group = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/tsv"
        os.makedirs(tsv_dir_group, exist_ok=True)

        for corr_param in corr_params:
            # ── average KDE density grids cell-by-cell across subjects ───────────────
            df_all_list = []
            for sub in subjects_to_group:
                tsv_dir_indiv = f"{main_dir}/{project_dir}/derivatives/pp_data/{sub}/{format_}/pmf/tsv"
                kde_fn = f"{tsv_dir_indiv}/{sub}_{cond_a}_vs_{cond_b}_{corr_param}-corr.tsv"
                if not os.path.isfile(kde_fn):
                    print(f"[SKIP] missing {kde_fn} — run make_corr_tsv.py for {sub} first")
                    continue
                df_all_list.append(pd.read_table(kde_fn))
            df_all = pd.concat(df_all_list, ignore_index=True)
            df_group = df_all.groupby(["roi", "x_grid", "y_grid"], sort=False)[["density"]].mean().reset_index()

            kde_fn_group = f"{tsv_dir_group}/{subject}_{cond_a}_vs_{cond_b}_{corr_param}-corr.tsv"
            df_group.to_csv(kde_fn_group, sep="\t", na_rep="NaN", index=False)
            print(f"Saving tsv: {kde_fn_group}")

            # ── average regression stats across subjects, with between-subject CI ────
            df_reg_all_list = []
            for sub in subjects_to_group:
                tsv_dir_indiv = f"{main_dir}/{project_dir}/derivatives/pp_data/{sub}/{format_}/pmf/tsv"
                reg_fn = f"{tsv_dir_indiv}/{sub}_{cond_a}_vs_{cond_b}_{corr_param}-regression.tsv"
                if not os.path.isfile(reg_fn):
                    print(f"[SKIP] missing {reg_fn} — run make_corr_tsv.py for {sub} first")
                    continue
                df_reg_all_list.append(pd.read_table(reg_fn))
            df_reg_all = pd.concat(df_reg_all_list, ignore_index=True)

            df_reg_group = df_reg_all.groupby(["roi"], sort=False)[
                ["slope", "intercept", "r2_pearson", "n_vertex", "mean_r2_weight"]
            ].mean().reset_index()
            for col in ["slope", "intercept", "r2_pearson"]:
                df_reg_group[f"{col}_ci_lower"] = df_reg_all.groupby(["roi"], sort=False)[col].quantile(0.025).values
                df_reg_group[f"{col}_ci_upper"] = df_reg_all.groupby(["roi"], sort=False)[col].quantile(0.975).values

            reg_fn_group = f"{tsv_dir_group}/{subject}_{cond_a}_vs_{cond_b}_{corr_param}-regression.tsv"
            df_reg_group.to_csv(reg_fn_group, sep="\t", na_rep="NaN", index=False)
            print(f"Saving tsv: {reg_fn_group}")

