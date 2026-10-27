"""
-----------------------------------------------------------------------------------------
check_medians.py
-----------------------------------------------------------------------------------------
Goal of the script:
Diagnostic only (writes nothing). Compare the leave-one-out median files already
saved (computed with the running pairwise median) with the true median across folds,
to know which outputs really change before re-running the pipeline.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: analysis name (e.g prf)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python check_medians.py /scratch/mszinte/data RetinoMaps prf
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import re
import sys
import glob
import numpy as np
from scipy import stats

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import load_surface
from maths_utils import multipletests_surface
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
analysis_name = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings = load_settings([os.path.join(base_dir, project_dir, "settings.yml"),
                          os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")])
analysis_info = settings[0]
subjects = analysis_info['subjects']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
task_name = analysis_info['analysis_task_names'][0]
dm_name = analysis_info['dm_name']
fdr_alpha = analysis_info['fdr_alpha']
maps_deriv = analysis_info['maps_names_css_loo']
pp_dir = f"{main_dir}/{project_dir}/derivatives/pp_data"

slope_idx, stderr_idx, trs_idx, rvalue_idx = 0, 4, 5, 2

def true_median(fns):
    return np.nanmedian(np.stack([load_surface(fn)[1] for fn in fns], axis=0), axis=0)

def sig_mask(stats_data):
    t = stats_data[slope_idx] / stats_data[stderr_idx]
    df = np.nanmax(stats_data[trs_idx]) - 2
    p = 2 * (1 - stats.t.cdf(np.abs(t), df=df))
    p_corr = multipletests_surface(p, correction='fdr_tsbh', alpha=fdr_alpha)
    return p_corr[0] <= 0.05

def n_sig(stats_data):
    return np.nansum(sig_mask(stats_data))

def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 2 else np.nan

for subject in subjects:
    for format_, extension in zip(formats, extensions):
        deriv_dir = f"{pp_dir}/{subject}/{format_}/{analysis_name}/prf_derivatives"
        hemis = ["_hemi-L", "_hemi-R"] if format_ == 'fsnative' else [""]
        for hemi in hemis:
            label = f"{subject} {format_}{hemi}"
            good = None

            # pRF stats
            fold_fns = sorted(glob.glob(f"{deriv_dir}/*task-{task_name}{hemi}_*loo-avg-*_{analysis_name}-css{dm_name}_stats.{extension}"))
            saved_fns = glob.glob(f"{deriv_dir}/*task-{task_name}{hemi}_*_loo-avg_{analysis_name}-css{dm_name}_stats.{extension}")
            if format_ != 'fsnative':
                fold_fns = [fn for fn in fold_fns if 'hemi-' not in fn]
            if len(fold_fns) > 0 and len(saved_fns) == 1:
                old = load_surface(saved_fns[0])[1]
                new = true_median(fold_fns)
                good = sig_mask(new)
                print(f"{label} | STATS  | folds={len(fold_fns)} | r(rvalue old,new)={corr(old[rvalue_idx], new[rvalue_idx]):.3f}"
                      f" | n sig FDR5% old={n_sig(old):.0f} new={n_sig(new):.0f}")

            # pRF derivatives
            fold_fns = sorted(glob.glob(f"{deriv_dir}/*task-{task_name}{hemi}_*loo-avg-*_{analysis_name}-css{dm_name}_deriv.{extension}"))
            saved_fns = glob.glob(f"{deriv_dir}/*task-{task_name}{hemi}_*_loo-avg_{analysis_name}-css{dm_name}_deriv.{extension}")
            if format_ != 'fsnative':
                fold_fns = [fn for fn in fold_fns if 'hemi-' not in fn]
            if len(fold_fns) > 0 and len(saved_fns) == 1:
                old = load_surface(saved_fns[0])[1]
                new = true_median(fold_fns)
                txt = []
                for map_name in ['prf_loo_rsq', 'prf_ecc', 'prf_size', 'prf_n']:
                    i = maps_deriv.index(map_name)
                    txt.append(f"{map_name}: r={corr(old[i], new[i]):.3f} med|diff|={np.nanmedian(np.abs(old[i]-new[i])):.3f}")
                print(f"{label} | DERIV  | folds={len(fold_fns)} | all vertices | " + " | ".join(txt))

                # Same comparison restricted to vertices kept in the analyses (new median, FDR 5%, ecc/size range)
                if good is not None:
                    ecc_i, size_i = maps_deriv.index('prf_ecc'), maps_deriv.index('prf_size')
                    keep = good & (new[ecc_i] >= analysis_info['ecc_th'][0]) & (new[ecc_i] <= analysis_info['ecc_th'][1]) \
                                & (new[size_i] >= analysis_info['size_th'][0]) & (new[size_i] <= analysis_info['size_th'][1])
                    txt = []
                    for map_name in ['prf_loo_rsq', 'prf_ecc', 'prf_size', 'prf_n']:
                        i = maps_deriv.index(map_name)
                        txt.append(f"{map_name}: r={corr(old[i][keep], new[i][keep]):.3f} med|diff|={np.nanmedian(np.abs(old[i][keep]-new[i][keep])):.3f}")
                    print(f"{label} | DERIV  | folds={len(fold_fns)} | kept n={keep.sum()} | " + " | ".join(txt))

        # GLM stats (2 folds: running median is correct, unless the median file was re-included or is stale)
        glm_dir = f"{pp_dir}/{subject}/{format_}/glm/glm_derivatives"
        saved_fns = sorted(glob.glob(f"{glm_dir}/*_loo-avg_glm-stats.{extension}"))
        for saved_fn in saved_fns:
            fold_fns = sorted(glob.glob(saved_fn.replace('_loo-avg_glm-stats', '_loo-avg-*_glm-stats')))
            name = os.path.basename(saved_fn)
            if len(fold_fns) != 2:
                print(f"{subject} {format_} | GLM    | {name} | WARNING: {len(fold_fns)} fold files found")
                continue
            old = load_surface(saved_fn)[1]
            new = true_median(fold_fns)
            same = np.allclose(old[[0, 1, 2, 4]], new[[0, 1, 2, 4]], equal_nan=True, rtol=1e-4, atol=1e-5)
            max_diff = np.nanmax(np.abs(old[rvalue_idx] - new[rvalue_idx]))
            print(f"{subject} {format_} | GLM    | {name} | identical: {same} | r(rvalue old,new)={corr(old[rvalue_idx], new[rvalue_idx]):.3f}"
                  f" | max|diff rvalue|={max_diff:.4f} | n sig FDR5% old={n_sig(old):.0f} new={n_sig(new):.0f}"
                  f" | median file newer than folds: {os.path.getmtime(saved_fn) >= max(os.path.getmtime(f) for f in fold_fns)}")
