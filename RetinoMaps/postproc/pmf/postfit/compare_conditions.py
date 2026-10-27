"""
-----------------------------------------------------------------------------------------
compare_conditions.py
-----------------------------------------------------------------------------------------
Goal of the script:
Mask a pRF derivatives image using the vertex mask produced by make_intertask_img.py,
combining a chosen set of intertask categories
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. pmf)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python compare_conditions.py /scratch/mszinte/data RetinoMaps pmf
-----------------------------------------------------------------------------------------
"""
import glob
import os
from itertools import combinations
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import bootstrap, permutation_test
from statsmodels.stats.multitest import multipletests

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import make_surface_image, load_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
analysis_name = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
task_names = analysis_info['analysis_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
averaging_templates = analysis_info['averaging_templates']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]
rois_method_format = "rois-group-mmp"
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
    "mPCS": "rgb(255, 111, 0)"
}

conditions = [
    ("pmf-css-rdm", "concat"),
    ("pmf-css-odm", "concat"),
    ("pmf-css-odm", "concat-residuals"),
]

early_rois = ["V1", "V2", "V3"]
higher_rois = ["iIPS", "sIPS", "sPCS", "iPCS", "mPCS"]


DM_ALPHA = {
    "pmf-css-rdm": 1.0,
    "pmf-css-odm": 0.55,
}
AVG_HATCH = {
    "concat": "",
    "concat-residuals": "//",
}

VALUE_COL = "prf_rsq_weighted_median"
N_RESAMPLES = 10000
FDR_ALPHA = 0.05

BAR_COLORS = [
    "rgba(165, 209, 216, 0.4)",
    "rgba(24, 151, 178, 0.4)",
]
X_POSITIONS = [0, 1]
BAR_WIDTH = 0.35

def rgb_to_mpl(rgb_str):
    r, g, b = map(int, rgb_str.replace("rgb(", "").replace(")", "").split(","))
    return (r / 255, g / 255, b / 255)


def condition_label(dm, am):
    return f"{dm}/{am}"

def summarize_across_subjects(df, value_col=VALUE_COL, confidence=0.95, n_resamples=N_RESAMPLES, seed=42):
    rows = []
    for (roi, dm, am), g in df.groupby(["roi", "dm", "avg_method"]):
        vals = g[value_col].dropna().values
        med = np.median(vals)
        if len(vals) > 1:
            res = bootstrap(
                (vals,), np.median,
                confidence_level=confidence,
                n_resamples=n_resamples,
                method="basic",
                random_state=seed,
            )
            ci_low, ci_high = res.confidence_interval.low, res.confidence_interval.high
        else:
            ci_low = ci_high = med
        rows.append({
            "roi": roi, "dm": dm, "avg_method": am,
            VALUE_COL: med,
            "prf_rsq_ci_down": ci_low,
            "prf_rsq_ci_up": ci_high,
            "n_subjects": len(vals),
        })
    return pd.DataFrame(rows)

def make_group_plot(data_sub, out_path):
    summary = summarize_across_subjects(data_sub)
    present_rois = set(data_sub["roi"].unique())
    rois = [roi for roi in roi_colors if roi in present_rois]

    x = np.arange(len(rois))
    n_cond = len(conditions)
    width = 0.8 / n_cond

    plt.figure(figsize=(13, 5))

    for i, (dm, am) in enumerate(conditions):
        sub = (
            summary[(summary["dm"] == dm) & (summary["avg_method"] == am)]
            .set_index("roi")
            .reindex(rois)
        )

        y = sub[VALUE_COL].values
        yerr = np.vstack([
            y - sub["prf_rsq_ci_down"].values,
            sub["prf_rsq_ci_up"].values - y,
        ])

        for j, roi in enumerate(rois):
            color = rgb_to_mpl(roi_colors[roi])
            plt.bar(
                x[j] + i * width,
                y[j],
                width,
                color=color,
                alpha=DM_ALPHA[dm],
                hatch=AVG_HATCH[am],
                edgecolor="black",
                yerr=[[yerr[0, j]], [yerr[1, j]]],
                capsize=3,
                linewidth=0.5,
            )

    plt.xticks(x + width * (n_cond - 1) / 2, rois, rotation=45, ha="right")
    plt.ylabel("Median R²")
    plt.title(
        f"Median R² per ROI × condition (task-{task_name}, n={data_sub['subject'].nunique()} subjects)"
        "\n(color = ROI, alpha = dm, hatch = avg_method; error bars = bootstrap 95% CI across subjects)"
    )

    legend_handles = [
        mpatches.Patch(
            facecolor="grey", alpha=DM_ALPHA[dm], hatch=AVG_HATCH[am], edgecolor="black",
            label=condition_label(dm, am),
        )
        for dm, am in conditions
    ]
    plt.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    print(f"Saved group plot -> {out_path}")
    return summary

def paired_permutation_test_between_conditions(df, value_col, cond_a, cond_b):
    """cond_a / cond_b are (dm, avg_method) tuples. Pairs observations by subject."""
    results = []
    dm_a, am_a = cond_a
    dm_b, am_b = cond_b
    for roi, group in df.groupby("roi"):
        a_rows = group[(group["dm"] == dm_a) & (group["avg_method"] == am_a)][["subject", value_col]]
        b_rows = group[(group["dm"] == dm_b) & (group["avg_method"] == am_b)][["subject", value_col]]
        merged = a_rows.merge(b_rows, on="subject", suffixes=("_a", "_b")).dropna()
        if merged.shape[0] < 2:
            continue
        a = merged[f"{value_col}_a"].values
        b = merged[f"{value_col}_b"].values
        res = permutation_test(
            (a, b),
            statistic=lambda x, y: (x - y).mean(),
            permutation_type="samples",
            alternative="two-sided",
            n_resamples=N_RESAMPLES,
        )
        results.append({
            "roi": roi,
            "condition_a": condition_label(*cond_a),
            "condition_b": condition_label(*cond_b),
            "mean_diff": (a - b).mean(),
            "p_value": res.pvalue,
            "n": len(merged),
        })
    return pd.DataFrame(results)





def plot_roi_comparison(data, roi_list, group_col, group_keys, group_labels, title, out_path, legend_title):
    pivots, stats = {}, {}
    for roi in roi_list:
        df_plot = data[(data["roi"] == roi) & (data[group_col].isin(group_keys))]
        df_pivot = (
            df_plot.pivot(index="subject", columns=group_col, values=VALUE_COL)[group_keys].dropna()
        )
        pivots[roi] = df_pivot
        stats[roi] = {k: (df_pivot[k].mean(), df_pivot[k].sem()) for k in group_keys}

    y_range = [0, 0.35]

    fig = make_subplots(rows=1, cols=len(roi_list), subplot_titles=[f"<b>{r}</b>" for r in roi_list])
    shown_subjects = set()

    for j, roi in enumerate(roi_list):
        col = j + 1
        df_pivot = pivots[roi]

        for k, (key, bar_color) in enumerate(zip(group_keys, BAR_COLORS)):
            mean_val, stderr_val = stats[roi][key]
            fig.add_trace(go.Bar(
                x=[X_POSITIONS[k]],
                y=[mean_val],
                error_y=dict(type="data", array=[stderr_val]),
                marker_color=bar_color,
                marker_line=dict(color=bar_color.replace("0.4", "0.9"), width=1.5),
                width=BAR_WIDTH,
                name=group_labels[k],
                showlegend=(j == 0),
                legendgroup=key,
            ), row=1, col=col)

        for subject, row_data in df_pivot.iterrows():
            show_legend = subject not in shown_subjects
            shown_subjects.add(subject)
            fig.add_trace(go.Scatter(
                x=X_POSITIONS,
                y=[row_data[k] for k in group_keys],
                mode="lines+markers",
                line=dict(color="rgba(100, 100, 100, 0.35)", width=1.5),
                marker=dict(size=8, color="rgba(100, 100, 100, 0.5)"),
                name=subject,
                showlegend=False,
                legendgroup=subject,
            ), row=1, col=col)

    fig.update_layout(
        height=550,
        width=max(350 * len(roi_list), 500),
        title_text=title,
        template="simple_white",
        showlegend=True,
        font=dict(size=16, family="Arial"),
        margin=dict(t=120),
        legend=dict(title=legend_title, orientation="v", x=1.02, y=0.5),
    )
    for col in range(1, len(roi_list) + 1):
        fig.update_yaxes(range=y_range, title_text="Median R²" if col == 1 else "", row=1, col=col)
        fig.update_xaxes(tickvals=X_POSITIONS, ticktext=group_labels, row=1, col=col)

    fig.write_image(out_path)
    print(f"Saved comparison plot -> {out_path}")
    return fig


def make_final_comparison_plots(data_sub):
    data_concat = data_sub[data_sub["avg_method"] == "concat"]
    model_keys = ["pmf-css-rdm", "pmf-css-odm"]
    model_labels = ["CSS RDM", "CSS ODM"]

    plot_roi_comparison(
        data_concat, early_rois, "dm", model_keys, model_labels,
        "CSS RDM vs CSS ODM (concat) — Early Visual ROIs",
        os.path.join(out_path, "concat_model_comparison_sacloc_early_rois.pdf"),
        legend_title="Model",
    )
    plot_roi_comparison(
        data_concat, higher_rois, "dm", model_keys, model_labels,
        "CSS RDM vs CSS ODM (concat) — Higher Visual ROIs",
        os.path.join(out_path, "concat_model_comparison_sacloc_higher_rois.pdf"),
        legend_title="Model",
    )

    data_odm = data_sub[data_sub["dm"] == "pmf-css-odm"]
    avg_method_keys = ["concat", "concat-residuals"]
    avg_method_labels = ["Concat", "Concat + Residuals"]

    plot_roi_comparison(
        data_odm, early_rois, "avg_method", avg_method_keys, avg_method_labels,
        "CSS ODM: Concat vs Concat + Residuals — Early Visual ROIs",
        os.path.join(out_path, "concat-residuals_model_comparison_sacloc_early_rois.pdf"),
        legend_title="Averaging method",
    )
    plot_roi_comparison(
        data_odm, higher_rois, "avg_method", avg_method_keys, avg_method_labels,
        "CSS ODM: Concat vs Concat + Residuals — Higher Visual ROIs",
        os.path.join(out_path, "concat-residuals_model_comparison_sacloc_higher_rois.pdf"),
        legend_title="Averaging method",
    )


# ── 1. load per-subject data for the 3 conditions ───────────────────────────
for format_ in formats: 
    task_name = task_names[0]
    out_path = f"{main_dir}/{project_dir}/derivatives/pp_data/group/{format_}/pmf/figures"
    subject_dirs = sorted(glob.glob(f"{main_dir}/{project_dir}/derivatives/pp_data/sub-*/{format_}"))
    print(subject_dirs)
    print(f"Found {len(subjects)} subjects")

    dfs = []
    for subject in subjects:
        print(f"loading data from subject {subject}")
        for dm, am in conditions:
            fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                task_name, preproc_prep, filtering, normalization, am, rois_method_format
            )
            path = (f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/tsv/"
                    f"{subject}_{fn_spec}_{dm}_params-median.tsv")
            if not os.path.exists(path):
                print(f"  Missing: {path}")
                continue
            df = pd.read_csv(path, sep="\t")
            df["subject"] = subject
            df["dm"] = dm
            df["avg_method"] = am
            dfs.append(df)

    data_sub = pd.concat(dfs, ignore_index=True)
    print(
        f"Loaded {len(data_sub)} rows | subjects: {data_sub['subject'].nunique()} "
        f"| rois: {list(data_sub['roi'].unique())}"
    )
# ── 2. group plot: median + bootstrap 95% CI per roi x condition ─────────── 
    make_group_plot(data_sub, os.path.join(out_path, "group_plot_sacloc.pdf"))
    
# ── 3. significance testing between the 3 conditions ────────────────────────
    all_results = []
    for cond_a, cond_b in combinations(conditions, 2):
        res = paired_permutation_test_between_conditions(data_sub, VALUE_COL, cond_a, cond_b)
        if res.empty:
            continue
        _, res["p_fdr"], _, _ = multipletests(res["p_value"], method="fdr_bh")
        res["significant"] = res["p_fdr"] < FDR_ALPHA
        all_results.append(res)

        print(f"\n── {condition_label(*cond_a)}  vs  {condition_label(*cond_b)} ──")
        print(
            res[["roi", "mean_diff", "p_value", "p_fdr", "significant", "n"]]
            .sort_values("p_fdr")
            .to_string(index=False)
        )

    results_df = pd.concat(all_results, ignore_index=True)
    #results_df.to_csv(out_path, index=False)
    print(f"\nSaved significance table -> {out_path}")

# ── 4. final comparison plots (bar + per-subject lines, one subplot per ROI) ─
    make_final_comparison_plots(data_sub)