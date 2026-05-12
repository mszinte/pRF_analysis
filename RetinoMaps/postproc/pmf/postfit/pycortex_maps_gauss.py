"""
-----------------------------------------------------------------------------------------
pycortex_maps_gauss.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset for gauss fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save in svg (e.g. no)
sys.argv[5]: fit type ('residuals', 'bold', or 'prf')
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON PERSONAL SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
>> python pycortex_maps_gauss.py [main directory] [project name] 
                                 [subject num] [save_in_svg] [fit_type]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python pycortex_maps_gauss.py ~/disks/meso_shared RetinoMaps sub-01 n residuals
python pycortex_maps_gauss.py ~/disks/meso_shared RetinoMaps sub-01 n bold
python pycortex_maps_gauss.py ~/disks/meso_shared RetinoMaps sub-01 n prf
python pycortex_maps_gauss.py ~/disks/meso_shared RetinoMaps hcp1.6mm n residuals
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
adapted by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import cortex
import importlib
import numpy as np
import matplotlib.pyplot as plt

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg = sys.argv[4]
fit_type = sys.argv[5]  # 'residuals', 'bold', or 'prf'

# Validate fit_type
valid_fit_types = ['residuals', 'bold', 'prf']
if fit_type not in valid_fit_types:
    sys.exit("Invalid fit_type '{}'. Must be one of: {}".format(fit_type, valid_fit_types))

# Define fit_type-specific filename patterns
fit_type_config = {
    'residuals': {
        'deriv_tag':  'pmf-residuals-gauss_deriv',
        'stats_tag':  'pmf-residuals-gauss_stats',
        'output_tag': 'gauss-residuals',
    },
    'bold': {
        'deriv_tag':  'pmf2-gauss_deriv',
        'stats_tag':  'pmf2-gauss_stats',
        'output_tag': 'gauss-bold',
    },
    'prf': {
        'deriv_tag':  'prf-gauss_deriv',
        'stats_tag':  'prf-gauss_stats',
        'output_tag': 'gauss-prf',
    },
}

cfg        = fit_type_config[fit_type]
deriv_tag  = cfg['deriv_tag']
stats_tag  = cfg['stats_tag']
output_tag = cfg['output_tag']

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_names = ['SacLoc']
maps_names_gauss = analysis_info['maps_names_gauss']
# Stats map names (shared with CSS pipeline: slope, intercept, rvalue, pvalue,
# stderr, trs, corr_pvalue_5pt, corr_pvalue_1pt)
maps_names_gauss_stats = analysis_info['maps_names_css_stats']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = ['concat']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# FDR alpha threshold
fdr_alpha = analysis_info['fdr_alpha']

# Plot scales
rsq_scale   = analysis_info['flatmap_rsq_scale']
ecc_scale   = analysis_info['flatmap_ecc_scale']
size_scale  = analysis_info['flatmap_size_scale']
alpha_range = analysis_info["flatmap_alpha_range"]

# Maps settings — build indices for deriv maps, then stats maps appended after
for idx, col_name in enumerate(maps_names_gauss):
    exec("{}_idx = idx".format(col_name))

# Stats indices are offset by the number of deriv maps
n_deriv_maps = len(maps_names_gauss)
for idx, col_name in enumerate(maps_names_gauss_stats):
    exec("{}_idx = {}".format(col_name, n_deriv_maps + idx))

cmap_polar, cmap_uni, cmap_ecc_size = 'hsv', 'Reds', 'Spectral'
col_offset = 1.0/14.0
cmap_steps = 255

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

for avg_method in avg_methods:
    if "loo" in avg_method:
        continue

    for format_ in formats:

        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:
            if 'drawn' in rois_method_format:
                overlay_fn = 'overlays.svg'
                rois_method_format_txt = ""
            else:
                overlay_fn = f"overlays_{rois_method_format}.svg"
                rois_method_format_txt = f"_{rois_method_format}"

            # Define directories
            prf_dir = "{}/{}/derivatives/pp_data/{}/{}/pmf".format(
                main_dir, project_dir, subject, format_)

            if not os.path.isdir(prf_dir):
                print(f"[SKIP] prf_dir not found for format={format_}: {prf_dir}")
                continue

            prf_deriv_dir  = "{}/pmf_derivatives".format(prf_dir)
            flatmaps_dir   = '{}/pycortex/flatmaps_{}'.format(prf_dir, output_tag)
            datasets_dir   = '{}/pycortex/datasets_{}'.format(prf_dir, output_tag)
            os.makedirs(flatmaps_dir, exist_ok=True)
            os.makedirs(datasets_dir, exist_ok=True)

            for prf_task_name in prf_task_names:

                # ------------------------------------------------------------------
                # Load derivative and stats data
                # ------------------------------------------------------------------
                if format_ == 'fsnative':
                    pycortex_subject = subject

                    deriv_avg_fn_L = '{}/{}_task-{}_hemi-L_{}_{}_{}_{}_{}.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, deriv_tag)
                    deriv_avg_fn_R = '{}/{}_task-{}_hemi-R_{}_{}_{}_{}_{}.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, deriv_tag)

                    stats_avg_fn_L = '{}/{}_task-{}_hemi-L_{}_{}_{}_{}_{}.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, stats_tag)
                    stats_avg_fn_R = '{}/{}_task-{}_hemi-R_{}_{}_{}_{}_{}.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, stats_tag)

                    print(f"Loading deriv: {os.path.basename(deriv_avg_fn_L)}")
                    deriv_mat = load_surface_pycortex(L_fn=deriv_avg_fn_L, R_fn=deriv_avg_fn_R)['data_concat']
                    print(f"Loading stats: {os.path.basename(stats_avg_fn_L)}")
                    stats_mat = load_surface_pycortex(L_fn=stats_avg_fn_L, R_fn=stats_avg_fn_R)['data_concat']

                elif format_ == '170k':
                    pycortex_subject = pycortex_subject_template

                    deriv_avg_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}.dtseries.nii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, deriv_tag)
                    stats_avg_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}.dtseries.nii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, stats_tag)

                    print(f"Loading deriv: {os.path.basename(deriv_avg_fn)}")
                    deriv_mat = load_surface_pycortex(brain_fn=deriv_avg_fn)['data_concat']
                    print(f"Loading stats: {os.path.basename(stats_avg_fn)}")
                    stats_mat = load_surface_pycortex(brain_fn=stats_avg_fn)['data_concat']

                # Combine deriv and stats into one matrix so indices align
                all_deriv_mat = np.concatenate((deriv_mat, stats_mat))

                # ------------------------------------------------------------------
                # Threshold data
                # ------------------------------------------------------------------
                amp_down     = all_deriv_mat[amplitude_idx,...]  > 0
                rsqr_th_down = all_deriv_mat[prf_rsq_idx,...]   >= analysis_info['rsqr_th']
                size_th_down = all_deriv_mat[prf_size_idx,...]  >= analysis_info['size_th'][0]
                size_th_up   = all_deriv_mat[prf_size_idx,...]  <= analysis_info['size_th'][1]
                ecc_th_down  = all_deriv_mat[prf_ecc_idx,...]   >= analysis_info['ecc_th'][0]
                ecc_th_up    = all_deriv_mat[prf_ecc_idx,...]   <= analysis_info['ecc_th'][1]

                # FDR-corrected p-value threshold (matches CSS script logic)
                if analysis_info['stats_th'] == 0.05:
                    stats_th_down = all_deriv_mat[corr_pvalue_5pt_idx,...] <= 0.05
                elif analysis_info['stats_th'] == 0.01:
                    stats_th_down = all_deriv_mat[corr_pvalue_1pt_idx,...] <= 0.01

                all_th = np.array((amp_down, rsqr_th_down, size_th_down,
                                   size_th_up, ecc_th_down, ecc_th_up, stats_th_down))
                all_deriv_mat[prf_rsq_idx, np.logical_and.reduce(all_th) == False] = 0

                # ------------------------------------------------------------------
                # Create flatmaps
                # ------------------------------------------------------------------
                print('Creating flatmaps ({})...'.format(output_tag))
                maps_names = []

                # r-square
                rsq_data = all_deriv_mat[prf_rsq_idx,...]
                alpha = (rsq_data - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
                alpha[alpha > 1] = 1
                param_rsq = {'data': rsq_data,
                             'cmap': cmap_uni,
                             'alpha': alpha,
                             'vmin': rsq_scale[0],
                             'vmax': rsq_scale[1],
                             'cbar': 'discrete',
                             'cortex_type': 'VertexRGB',
                             'description': 'Gaussian pRF R2',
                             'curv_brightness': 1,
                             'curv_contrast': 0.1,
                             'add_roi': save_svg,
                             'cbar_label': 'pRF R2',
                             'overlay_fn': overlay_fn,
                             'with_labels': True}
                maps_names.append('rsq')

                # polar angle
                pol_comp_num = all_deriv_mat[polar_real_idx,...] + 1j * all_deriv_mat[polar_imag_idx,...]
                polar_ang = np.angle(pol_comp_num)
                ang_norm  = (polar_ang + np.pi) / (np.pi * 2.0)
                ang_norm  = np.fmod(ang_norm + col_offset, 1)
                param_polar = {'data': ang_norm,
                               'cmap': cmap_polar,
                               'alpha': alpha,
                               'vmin': 0,
                               'vmax': 1,
                               'cmap_steps': cmap_steps,
                               'cortex_type': 'VertexRGB',
                               'cbar': 'polar',
                               'col_offset': col_offset,
                               'description': 'Gaussian pRF polar angle',
                               'curv_brightness': 0.1,
                               'curv_contrast': 0.25,
                               'add_roi': save_svg,
                               'overlay_fn': overlay_fn,
                               'with_labels': True}
                exec('param_polar_{} = param_polar'.format(int(cmap_steps)))
                exec('maps_names.append("polar_{}")'.format(int(cmap_steps)))

                # eccentricity
                ecc_data = all_deriv_mat[prf_ecc_idx,...]
                param_ecc = {'data': ecc_data,
                             'cmap': cmap_ecc_size,
                             'alpha': alpha,
                             'vmin': ecc_scale[0],
                             'vmax': ecc_scale[1],
                             'cbar': 'ecc',
                             'cortex_type': 'VertexRGB',
                             'description': 'Gaussian pRF eccentricity',
                             'curv_brightness': 1,
                             'curv_contrast': 0.1,
                             'add_roi': save_svg,
                             'overlay_fn': overlay_fn,
                             'with_labels': True}
                maps_names.append('ecc')

                # size
                size_data = all_deriv_mat[prf_size_idx,...]
                param_size = {'data': size_data,
                              'cmap': cmap_ecc_size,
                              'alpha': alpha,
                              'vmin': size_scale[0],
                              'vmax': size_scale[1],
                              'cbar': 'discrete',
                              'cortex_type': 'VertexRGB',
                              'description': 'Gaussian pRF size',
                              'curv_brightness': 1,
                              'curv_contrast': 0.1,
                              'add_roi': False,
                              'cbar_label': 'pRF size',
                              'overlay_fn': overlay_fn,
                              'with_labels': True}
                maps_names.append('size')

                # ------------------------------------------------------------------
                # Draw flatmaps and build dataset
                # ------------------------------------------------------------------
                volumes = {}
                for maps_name in maps_names:

                    roi_name  = 'prf_{}'.format(maps_name)
                    roi_param = {'subject': pycortex_subject, 'roi_name': roi_name}
                    print(roi_name)
                    exec('param_{}.update(roi_param)'.format(maps_name))
                    exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(
                        maps_name=maps_name))

                    # PDF filename uses output_tag
                    exec("plt.savefig('{}/{}_task-{}_{}_{}_{}_{}{}_{}_{}.pdf')".format(
                        flatmaps_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method,
                        rois_method_format_txt, output_tag, maps_name))
                    plt.close()

                    exec('vol_description = param_{}["description"]'.format(maps_name))
                    exec('volume = volume_{}'.format(maps_name))
                    volumes.update({vol_description: volume})

                # HDF dataset filename uses output_tag
                dataset_file = "{}/{}_task-{}_{}_{}_{}_{}{}_{}_{}.hdf".format(
                    datasets_dir, subject, prf_task_name,
                    preproc_prep, filtering, normalization,
                    avg_method, rois_method_format_txt, output_tag, fit_type)
                if os.path.exists(dataset_file):
                    os.system("rm -fv {}".format(dataset_file))
                dataset = cortex.Dataset(data=volumes)
                dataset.save(dataset_file)