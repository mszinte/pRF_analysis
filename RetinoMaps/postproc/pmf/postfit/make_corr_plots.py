"""
-----------------------------------------------------------------------------------------
make_corr_plots.py
-----------------------------------------------------------------------------------------
Goal of the script:
Plot the ROI x parameter grid of KDE-density correlation panels (position, size, R2)
between two SacLoc / pRF conditions, for either a single subject or a group of subjects.
Reads the TSVs produced by make_corr_tsv.py — this script does not touch any raw surface
data, so it runs identically for an individual subject or a group (a group's TSVs were
already averaged across subjects by make_corr_tsv.py).

Rows = ROIs, columns = parameters. x-axis = condition A, y-axis = condition B. Each panel
shows: a dashed identity line (x = y) and a filled KDE contour of the (weighted) joint
density of vertex-level parameter estimates for that ROI.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: subject (e.g. sub-11) OR group name (e.g. group) — must match whatever
             was passed to make_corr_tsv.py, whose TSVs this script reads
sys.argv[4]: condition A (e.g. pmf-css-rdm)
sys.argv[5]: condition B (e.g. pmf-css-odm, or prf-css for a pmf-vs-prf comparison)
-----------------------------------------------------------------------------------------
Output(s):
{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format}/pmf/figures/
    {subject}_{cond_a}_vs_{cond_b}_corr.pdf
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python make_corr_plots.py /scratch/mszinte/data RetinoMaps sub-11 pmf-css-odm prf-css
python make_corr_plots.py /scratch/mszinte/data RetinoMaps group pmf-css-odm prf-css
-----------------------------------------------------------------------------------------
Adapted from Martin Szinte's make_corr_fig.py (amblyo7T_prf)
-----------------------------------------------------------------------------------------
"""
import os
import sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
cond_a = sys.argv[4]
cond_b = sys.argv[5]

# Must match the ROI names / corr_params used in make_corr_tsv.py
roi_colors = {
    "V1": "rgb(243, 231, 155)",
    "V2": "rgb(250, 196, 132)",
    "V3": "rgb(248, 160, 126)",
    "V3AB": "rgb(235, 127, 134)",
    "LO": "rgb(150, 0, 90)",
    "VO": "rgb(0, 0, 200)",
    "hMT+": "rgb(0, 25, 255)",
    "iIPS": "rgb(0, 152, 255)",
    "sIPS": "rgb(44, 255, 150)",
    "iPCS": "rgb(151, 255, 0)",
    "sPCS": "rgb(255, 234, 0)",
    "mPCS": "rgb(255, 111, 0)",
}
rois = list(roi_colors.keys())

axis_labels = {
    "prf_x": "X coord. (dva)",
    "prf_y": "Y coord. (dva)",
    "prf_ecc": "ecc. (dva)",
    "prf_size": "size (dva)",
    "prf_rsq": "R²",
}

analysis_type_a, analysis_type_b = cond_a.split("-")[0], cond_b.split("-")[0]
label_a = cond_a.split("-")[-1].upper() if analysis_type_a == "pmf" else "pRF"
label_b = cond_b.split("-")[-1].upper() if analysis_type_b == "pmf" else "pRF"


base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_type_a}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
analysis_info = settings[0]
formats = analysis_info["formats"]
corr_params = analysis_info["corr_params"]
corr_param_ranges = {
    'prf_rsq':    analysis_info['prf_rsq_param_range'],
    'prf_x':      analysis_info['prf_x_param_range'],
    'prf_y':      analysis_info['prf_y_param_range'],
    'prf_size':   analysis_info['prf_size_param_range'],
    'prf_ecc':    analysis_info['prf_ecc_param_range']
}

# Figure layout settings
roi_colors = analysis_info['roi_colors']
fig_margin = analysis_info['rois_fig_margin']
#rois = analysis_info['rois']
rois_hor_spacing = analysis_info['rois_hor_spacing']
rois_ver_spacing = analysis_info['rois_ver_spacing']
rois_plot_height = analysis_info['rois_plot_height']
rois_plot_width = analysis_info['rois_plot_width']
rows, cols = len(rois), len(corr_params)

fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows - 1))
fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols - 1))
hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

cell_size = 320
margin = dict(l=90, r=20, t=20, b=90)


for format_ in formats:
    print(f"\n=== format: {format_} ===")
    tsv_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/tsv"
    fig_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/figures"
    if not os.path.isdir(tsv_dir):
        print(f"[SKIP] tsv_dir not found: {tsv_dir}")
        continue
    os.makedirs(fig_dir, exist_ok=True)

    # ── 1. load the KDE + regression tables for every parameter ─────────────────────
    df_kde_by_param = {}
    df_reg_parts = []
    for corr_param in corr_params:
        kde_fn = f"{tsv_dir}/{subject}_{cond_a}_vs_{cond_b}_{corr_param}-corr.tsv"
        reg_fn = f"{tsv_dir}/{subject}_{cond_a}_vs_{cond_b}_{corr_param}-regression.tsv"
        df_kde_by_param[corr_param] = pd.read_table(kde_fn)
        df_reg_p = pd.read_table(reg_fn)
        df_reg_p.insert(0, "corr_param", corr_param)
        df_reg_parts.append(df_reg_p)
    df_reg = pd.concat(df_reg_parts, ignore_index=True)

    print("\n── regression stats (not plotted, for reference) ──")
    print_cols = [c for c in ["corr_param", "roi", "slope", "r2_pearson", "n_vertex",
                               "slope_ci_lower", "slope_ci_upper"] if c in df_reg.columns]
    print(df_reg[print_cols].sort_values(["corr_param", "roi"]).to_string(index=False))

    # ── 2. build the ROI x parameter grid of KDE-density panels ─────────────────────
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, horizontal_spacing=hor_spacing, vertical_spacing=ver_spacing)

    for l, corr_param in enumerate(corr_params):
        df_kde_p = df_kde_by_param[corr_param]
        x_grid = np.sort(df_kde_p["x_grid"].unique())
        y_grid = np.sort(df_kde_p["y_grid"].unique())
        param_range = corr_param_ranges[corr_param]

        for j, roi in enumerate(rois):
            roi_color = roi_colors[roi]
            r, g, b = [int(v) for v in roi_color.replace("rgb(", "").replace(")", "").split(",")]
            df_roi = df_kde_p[df_kde_p["roi"] == roi].sort_values(["y_grid", "x_grid"])
            density = df_roi["density"].values.reshape(len(y_grid), len(x_grid))
            d_max = np.nanmax(density)
            density_norm = density / d_max if d_max > 0 else density

            line_xy = np.linspace(param_range[0], param_range[1], 50)
            fig.add_trace(go.Scatter(x=line_xy, y=line_xy, mode="lines",
                          line=dict(color="black", width=2, dash="dash"),
                          opacity=0.5, showlegend=False), row=j + 1, col=l + 1)
            fig.add_trace(go.Contour(
                x=x_grid, y=y_grid, z=density_norm,
                colorscale=[[0, f"rgba({r},{g},{b},0)"], [1, f"rgba({r},{g},{b},1)"]],
                showscale=False,
                contours=dict(coloring="fill", showlines=False, start=0.05, end=1.0, size=0.05),
                zmin=0, zmax=1), row=j + 1, col=l + 1)

            reg_row = df_reg[(df_reg["corr_param"] == corr_param) & (df_reg["roi"] == roi)]
            n_vertex = int(reg_row["n_vertex"].values[0]) if len(reg_row) and not np.isnan(reg_row["n_vertex"].values[0]) else 0
            fig.add_annotation(x=param_range[1] - 0.05 * (param_range[1] - param_range[0]),y=param_range[0] + 0.05 * (param_range[1] - param_range[0]),text=f"{roi}",
                               xanchor="right", yanchor="bottom",showarrow=False, font=dict(color=roi_color, size=14), row=j + 1, col=l + 1)
            x_title = f"{analysis_type_a} – {axis_labels[corr_param]}" if j == rows - 1 else ""
            fig.update_xaxes(title_text=x_title, range=param_range, showline=True, title_font=dict(color="black", family="Rockwell", size=20), row=j + 1, col=l + 1)
            fig.update_yaxes(title_text=f"{analysis_type_b} – {axis_labels[corr_param]}", range=param_range,showline=True, title_font=dict(color="black", family="Rockwell", size=20), title_standoff=0, row=j + 1, col=l + 1)
            idx = j * cols + l + 1
            x_ref = "x" if idx == 1 else f"x{idx}"
            fig.update_yaxes(scaleanchor=x_ref, scaleratio=1, row=j + 1, col=l + 1)

    fig.update_layout(
        template = "simple_white",
        width=fig_width + margin["l"] + margin["r"],
        height=fig_height + margin["t"] + margin["b"],
        margin=margin, showlegend=False, plot_bgcolor="white",font=dict(family="Rockwell", size=22),
    )

    out_path = os.path.join(fig_dir, f"{subject}_{cond_a}_vs_{cond_b}_corr.pdf")
    fig.write_image(out_path)
    print(f"Saved correlation plot -> {out_path}")
