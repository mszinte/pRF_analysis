"""
-----------------------------------------------------------------------------------------
make_fit_comparison_plots.py
-----------------------------------------------------------------------------------------
Goal of the script:
Plot the ROI x parameter grid of KDE-density comparison panels (position, size, n, R2)
between two CSS fits of the same data (e.g. prfpy vs prfmodel), for every fit found for
one subject (all formats x averaging methods x tasks in the settings).

Rows = ROIs, columns = parameters. x-axis = fit A, y-axis = fit B. Each panel shows a
dashed identity line (x = y) and a filled contour of the R2-weighted 2D KDE of the
vertex-level estimates for that ROI (vertices with R2 > RSQ_THRESHOLD in both fits,
weights = mean R2 of both fits, fixed GRID_SIZE x GRID_SIZE grid, density normalised
per panel), as in make_corr_tsv.py / make_corr_plots.py (RetinoMaps, pmf).
Also saves per ROI and parameter: n_vertex, Pearson r and r2, p-value, weighted Deming
slope and intercept, median difference (B - A) and mean R2 weight.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: subject name (e.g. sub-03)
sys.argv[4]: analysis name (e.g. prf)
sys.argv[5]: fit folder A (e.g. fit_prfanalysis, the prfpy fit)
sys.argv[6]: fit folder B (e.g. fit, the prfmodel fit)
-----------------------------------------------------------------------------------------
Output(s):
{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format}/{output_folder}/
    figures/{subject}_task-{task}_{avg_method}_{fit A}_vs_{fit B}_corr.pdf
    tsv/{subject}_task-{task}_{avg_method}_{fit A}_vs_{fit B}_corr.tsv
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_fit_comparison_plots.py [main directory] [project name] [subject]
                                    [analysis name] [fit folder A] [fit folder B]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_fit_comparison_plots.py /scratch/mszinte/data amsterdam24 sub-03 prf fit_prfanalysis fit
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
Adapted from make_corr_tsv.py / make_corr_plots.py (RetinoMaps)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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
analysis_name = sys.argv[4]
fit_folder_a = sys.argv[5]
fit_folder_b = sys.argv[6]

RSQ_THRESHOLD = 0.1
GRID_SIZE = 100
label_a, label_b = fit_folder_a, fit_folder_b

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
task_names = analysis_info['analysis_task_names']
avg_methods = analysis_info['avg_methods']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
output_folder = analysis_info['output_folder']
dm_name = analysis_info['dm_name']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']
max_ecc_size = analysis_info['max_ecc_size']
roi_colors = {roi: np.array(rgb) / 255 for roi, rgb in analysis_info['rois_colors'].items()}

# Parameters to compare: fit map name -> (axis label, axis range)
corr_params = {
    "mu_x":      ("X coord. (dva)", [-max_ecc_size, max_ecc_size]),
    "mu_y":      ("Y coord. (dva)", [-max_ecc_size, max_ecc_size]),
    "prf_ecc":   ("ecc. (dva)", [0, max_ecc_size]),
    "prf_size":  ("size (dva)", [0, max_ecc_size]),
    "n":         ("n", [0, 1.5]),
    "r_squared": ("R²", [0, 0.8]),
}

# Map names of the fit files written by prf_cssfit.py / prf_cssfit_gpu.py
maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
              'n', 'hrf_1', 'hrf_2', 'r_squared']

# Set pycortex db
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)


# ── helpers ─────────────────────────────────────────────────────────────────────────
def fit_fns(fit_dir, format_, extension, task_name, avg_method):
    """Fit file names as written by the fit scripts (fsnative: one file per hemisphere)."""
    base = f"{subject}_task-{task_name}_{{hemi}}{preproc_prep}_{filtering}_{normalization}_{avg_method}"
    suffix = f"{analysis_name}-css{dm_name}_fit.{extension}"
    if format_ == 'fsnative':
        return [f"{fit_dir}/{base.format(hemi=f'{hemi}_')}_{suffix}" for hemi in ['hemi-L', 'hemi-R']]
    return [f"{fit_dir}/{base.format(hemi='')}_{suffix}"]


def load_fit(fns):
    """Fit maps of all files concatenated over vertices -> DataFrame (one row per vertex)."""
    fit_mat = np.concatenate([load_surface(fn)[1] for fn in fns], axis=1)
    names = maps_names + ['loo_r_squared'] if fit_mat.shape[0] == len(maps_names) + 1 else maps_names
    params = pd.DataFrame(fit_mat.T, columns=names)
    params['prf_ecc'] = np.hypot(params['mu_x'], params['mu_y'])
    return params


def roi_labels(format_, rois_method, n_vertices):
    """ROI name per vertex, same ROI masks as make_tsv_css.py ('none' outside ROIs)."""
    pycortex_subject = subject if format_ == 'fsnative' else pycortex_subject_template
    labels = np.full(n_vertices, 'none', dtype=object)
    start = 0
    for hemi in ['hemi-L', 'hemi-R']:
        roi_verts = get_rois(subject=pycortex_subject, surf_format=format_, rois_type=rois_method,
                             mask=True, hemis=hemi)
        if format_ == 'fsnative':
            # masks are per hemisphere: left vertices first, then right
            n_hemi = len(next(iter(roi_verts.values())))
            for roi, roi_mask in roi_verts.items():
                labels[start:start + n_hemi][roi_mask] = roi
            start += n_hemi
        else:
            # masks index the whole 170k brain
            for roi, roi_mask in roi_verts.items():
                labels[roi_mask] = roi
    return labels


def weighted_kde(x, y, weights, param_range):
    """R2-weighted 2D Gaussian KDE on a fixed grid, normalised to a maximum of 1."""
    grid = np.linspace(param_range[0], param_range[1], GRID_SIZE)
    xx, yy = np.meshgrid(grid, grid)
    try:
        kde = gaussian_kde(np.vstack([x, y]), weights=weights / weights.sum())
    except np.linalg.LinAlgError:
        # e.g. all values identical in both fits (parameter at a bound): no 2D density
        return grid, None
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(GRID_SIZE, GRID_SIZE)
    return grid, density / density.max()


def weighted_deming_regression(x, y, weights=None):
    """Weighted Deming (orthogonal) regression, symmetric in x and y (from make_corr_tsv.py)."""
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


# ── loop over all fits found ────────────────────────────────────────────────────────
for format_, extension in zip(formats, extensions):
    prf_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{output_folder}"
    if not os.path.isdir(prf_dir):
        print(f"[SKIP] prf_dir not found for format={format_}: {prf_dir}")
        continue
    fig_dir, tsv_dir = f"{prf_dir}/figures", f"{prf_dir}/tsv"
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(tsv_dir, exist_ok=True)

    rois_method = rois_methods[format_][0]
    rois = (list(analysis_info[rois_method].keys()) if rois_method == 'rois-group-mmp'
            else analysis_info[rois_method])

    for avg_method in avg_methods:
        for task_name in task_names:
            print(f"\n=== {format_} - {avg_method} - {task_name} ===")

            # both fits must exist
            fns_a = fit_fns(f"{prf_dir}/{fit_folder_a}", format_, extension, task_name, avg_method)
            fns_b = fit_fns(f"{prf_dir}/{fit_folder_b}", format_, extension, task_name, avg_method)
            missing = [fn for fn in fns_a + fns_b if not os.path.isfile(fn)]
            if missing:
                print(f"[SKIP] not found in {prf_dir}: {[fn.replace(prf_dir + '/', '') for fn in missing]}")
                continue

            fit_a, fit_b = load_fit(fns_a), load_fit(fns_b)
            if len(fit_a) != len(fit_b):
                print(f"[SKIP] different number of vertices: {len(fit_a)} vs {len(fit_b)}")
                continue
            roi = roi_labels(format_, rois_method, len(fit_a))
            rois_found = [r for r in rois if r in set(roi)]

            rsq_a, rsq_b = fit_a['r_squared'].to_numpy(), fit_b['r_squared'].to_numpy()
            r2_combined = (rsq_a + rsq_b) / 2

            # ROI x parameter grid of KDE panels
            fig, axes = plt.subplots(len(rois_found), len(corr_params), squeeze=False,
                                     figsize=(2.6 * len(corr_params), 2.6 * len(rois_found)))
            stats_rows = []
            for j, roi_name in enumerate(rois_found):
                roi_color = roi_colors.get(roi_name, np.array([0.3, 0.3, 0.3]))
                for l, (param, (axis_label, param_range)) in enumerate(corr_params.items()):
                    ax = axes[j, l]
                    x_vals = fit_a[param].to_numpy()
                    y_vals = fit_b[param].to_numpy()

                    # ROI vertices, no NaN, R2 above threshold in both fits (KDE weights must be > 0)
                    mask = ((roi == roi_name) & ~(np.isnan(x_vals) | np.isnan(y_vals) | np.isnan(r2_combined))
                            & (rsq_a > RSQ_THRESHOLD) & (rsq_b > RSQ_THRESHOLD))
                    x_vals, y_vals, weights = x_vals[mask], y_vals[mask], r2_combined[mask]
                    n_vertex = len(x_vals)

                    # KDE panel
                    ax.plot(param_range, param_range, color="k", ls="--", lw=1, alpha=0.5)
                    if n_vertex > 2:
                        grid, density = weighted_kde(x_vals, y_vals, weights, param_range)
                        if density is not None:
                            cmap = LinearSegmentedColormap.from_list(roi_name, [(*roi_color, 0), (*roi_color, 1)])
                            ax.contourf(grid, grid, density, levels=np.arange(0.05, 1.0001, 0.05),
                                        cmap=cmap, vmin=0, vmax=1)
                    ax.set_xlim(param_range)
                    ax.set_ylim(param_range)
                    ax.set_aspect("equal")
                    ax.spines[["top", "right"]].set_visible(False)
                    ax.text(0.95, 0.05, roi_name, transform=ax.transAxes, ha="right", va="bottom",
                            color=roi_color, fontsize=10, fontweight="bold")
                    ax.set_ylabel(f"{label_b} – {axis_label}", fontsize=8)
                    if j == len(rois_found) - 1:
                        ax.set_xlabel(f"{label_a} – {axis_label}", fontsize=8)

                    # statistics
                    if n_vertex > 2:
                        slope, intercept = weighted_deming_regression(x_vals, y_vals, weights)
                        r_val, p_val = pearsonr(x_vals, y_vals)
                    else:
                        slope = intercept = r_val = p_val = np.nan
                    stats_rows.append({
                        "roi": roi_name, "corr_param": param, "n_vertex": n_vertex,
                        "r_pearson": r_val, "r2_pearson": r_val ** 2, "p_value": p_val,
                        "slope": slope, "intercept": intercept,
                        "median_diff": np.nanmedian(y_vals - x_vals) if n_vertex else np.nan,
                        "mean_r2_weight": np.nanmean(weights) if n_vertex else np.nan})

            fig.suptitle(f"{subject} · {format_} · {avg_method} · task-{task_name}: "
                         f"{label_a} (x) vs {label_b} (y), R² > {RSQ_THRESHOLD} in both",
                         x=0.01, ha="left", fontsize=12)
            plt.tight_layout(rect=(0, 0, 1, 0.98))

            fn_base = f"{subject}_task-{task_name}_{avg_method}_{label_a}_vs_{label_b}_corr"
            fig.savefig(f"{fig_dir}/{fn_base}.pdf")
            plt.close(fig)
            print(f"Saved figure: {fig_dir}/{fn_base}.pdf")

            df_stats = pd.DataFrame(stats_rows)
            df_stats.to_csv(f"{tsv_dir}/{fn_base}.tsv", sep="\t", na_rep="NaN", index=False)
            print(f"Saved tsv: {tsv_dir}/{fn_base}.tsv")
            print(df_stats.pivot(index="roi", columns="corr_param", values="r_pearson")
                  .loc[rois_found, list(corr_params)].round(2).to_string())
