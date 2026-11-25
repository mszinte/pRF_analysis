"""
-----------------------------------------------------------------------------------------
pycortex_maps_od.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset of the ocular dominance (OD) results
(output of make_od_img.py). Vertices are thresholded with the pRF analysis
thresholds applied to both eyes, and the transparency follows the mean pRF R2
of both eyes. Colour: left eye stronger = red (left of the colour bar),
right eye stronger = blue (right of the colour bar).
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. prf)
sys.argv[5]: save in svg (e.g. no)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_od.py [main directory] [project] [subject] [analysis name] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python pycortex_maps_od.py ~/disks/meso_shared amblyo7T_prf sub-01 prf n
python pycortex_maps_od.py ~/disks/meso_shared amblyo7T_prf sub-20 prf n
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import numpy as np
import matplotlib.pyplot as plt

# Personal import
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
save_svg = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
maps_names_css = analysis_info['maps_names_css']
maps_names_css_stats = analysis_info['maps_names_css_stats']
maps_names_css_od = analysis_info['maps_names_css_od']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalizations = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]
od_task_names = analysis_info['od_task_names']
eyes = ['LeftEye', 'RightEye']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Maps settings (applied to right - left: left eye = low values = red, right eye = high values = blue)
cmap_od = 'coolwarm_r'

# plot scales
od_scale = analysis_info['flatmap_od_scale']
od_z_scale = analysis_info['flatmap_od_z_scale']
alpha_range = analysis_info["flatmap_alpha_range"]

# Index of the maps
for idx, col_name in enumerate(maps_names_css + maps_names_css_stats):
    exec("{}_idx = idx".format(col_name))
for idx, col_name in enumerate(maps_names_css_od):
    exec("{}_idx = idx".format(col_name))

for normalization in normalizations:
    for avg_method in avg_methods:
        for format_ in formats:

            # define list of rois for each format
            rois_methods_format = rois_methods[format_]
            for rois_method_format in rois_methods_format:

                # Define directories and fn
                prf_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{output_folder}"

                if not os.path.isdir(prf_dir):
                    print(f"[SKIP] prf_dir not found for format={format_}: {prf_dir}")
                    continue

                prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
                flatmaps_dir = '{}/pycortex/flatmaps_od'.format(prf_dir)
                datasets_dir = '{}/pycortex/datasets_od'.format(prf_dir)
                os.makedirs(flatmaps_dir, exist_ok=True)
                os.makedirs(datasets_dir, exist_ok=True)

                # Define overlay
                overlay_fn = f"overlays_{rois_method_format}.svg"

                for task_name in od_task_names:

                    # Define files: OD image, derivatives and stats of each eye
                    od_fn_spec = f'{subject}_task-{task_name}{{hemi}}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}'
                    eye_fn_spec = f'{subject}_task-{task_name}{{eye}}{{hemi}}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}'

                    if format_ == 'fsnative':
                        pycortex_subject = subject
                        extension = 'func.gii'
                        od_fn_L = f"{prf_deriv_dir}/{od_fn_spec.format(hemi='_hemi-L')}_od.{extension}"
                        od_fn_R = f"{prf_deriv_dir}/{od_fn_spec.format(hemi='_hemi-R')}_od.{extension}"
                        if not os.path.isfile(od_fn_L) or not os.path.isfile(od_fn_R):
                            print(f"[SKIP] OD image not found: {od_fn_L}")
                            continue
                        print(f"Loading OD: {os.path.basename(od_fn_L)}")
                        od_mat = load_surface_pycortex(L_fn=od_fn_L, R_fn=od_fn_R)['data_concat']

                        deriv_mats, stats_mats = {}, {}
                        for eye in eyes:
                            deriv_fn_L = f"{prf_deriv_dir}/{eye_fn_spec.format(eye=eye, hemi='_hemi-L')}_deriv.{extension}"
                            deriv_fn_R = f"{prf_deriv_dir}/{eye_fn_spec.format(eye=eye, hemi='_hemi-R')}_deriv.{extension}"
                            stats_fn_L = deriv_fn_L.replace(f'{analysis_name}-css{dm_name}_deriv', f'{analysis_name}-css{dm_name}_stats')
                            stats_fn_R = deriv_fn_R.replace(f'{analysis_name}-css{dm_name}_deriv', f'{analysis_name}-css{dm_name}_stats')
                            print(f"Loading deriv: {os.path.basename(deriv_fn_L)}")
                            deriv_mats[eye] = load_surface_pycortex(L_fn=deriv_fn_L, R_fn=deriv_fn_R)['data_concat']
                            stats_mats[eye] = load_surface_pycortex(L_fn=stats_fn_L, R_fn=stats_fn_R)['data_concat']

                    elif format_ == '170k':
                        pycortex_subject = pycortex_subject_template
                        extension = 'dtseries.nii'
                        if subject != 'sub-170k': save_svg = False
                        od_fn = f"{prf_deriv_dir}/{od_fn_spec.format(hemi='')}_od.{extension}"
                        if not os.path.isfile(od_fn):
                            print(f"[SKIP] OD image not found: {od_fn}")
                            continue
                        print(f"Loading OD: {os.path.basename(od_fn)}")
                        od_mat = load_surface_pycortex(brain_fn=od_fn)['data_concat']

                        deriv_mats, stats_mats = {}, {}
                        for eye in eyes:
                            deriv_fn = f"{prf_deriv_dir}/{eye_fn_spec.format(eye=eye, hemi='')}_deriv.{extension}"
                            stats_fn = deriv_fn.replace(f'{analysis_name}-css{dm_name}_deriv', f'{analysis_name}-css{dm_name}_stats')
                            print(f"Loading deriv: {os.path.basename(deriv_fn)}")
                            deriv_mats[eye] = load_surface_pycortex(brain_fn=deriv_fn)['data_concat']
                            stats_mats[eye] = load_surface_pycortex(brain_fn=stats_fn)['data_concat']

                    # Threshold mat: pRF analysis thresholds applied to both eyes
                    all_th_eyes = []
                    for eye in eyes:
                        all_deriv_mat_th = np.concatenate([deriv_mats[eye], stats_mats[eye]])
                        amp_down = all_deriv_mat_th[amplitude_idx,...] > 0
                        rsq_down = all_deriv_mat_th[prf_rsq_idx,...] >= analysis_info['rsqr_th']
                        size_th_down = all_deriv_mat_th[prf_size_idx,...] >= analysis_info['size_th'][0]
                        size_th_up = all_deriv_mat_th[prf_size_idx,...] <= analysis_info['size_th'][1]
                        ecc_th_down = all_deriv_mat_th[prf_ecc_idx,...] >= analysis_info['ecc_th'][0]
                        ecc_th_up = all_deriv_mat_th[prf_ecc_idx,...] <= analysis_info['ecc_th'][1]
                        n_th_down = all_deriv_mat_th[prf_n_idx,...] >= analysis_info['n_th'][0]
                        n_th_up = all_deriv_mat_th[prf_n_idx,...] <= analysis_info['n_th'][1]
                        if analysis_info['stats_th'] == 0.05: stats_th_down = all_deriv_mat_th[corr_pvalue_5pt_idx,...] <= 0.05
                        elif analysis_info['stats_th'] == 0.01: stats_th_down = all_deriv_mat_th[corr_pvalue_1pt_idx,...] <= 0.01
                        all_th_eyes.append(np.logical_and.reduce(np.array((amp_down, rsq_down, size_th_down, size_th_up,
                                                                           ecc_th_down, ecc_th_up, n_th_down, n_th_up,
                                                                           stats_th_down))))
                    all_th = np.logical_and(all_th_eyes[0], all_th_eyes[1])

                    # Alpha: mean R2 of both eyes, 0 outside the thresholds
                    rsq_data = np.nanmean([deriv_mats[eye][prf_rsq_idx,...] for eye in eyes], axis=0)
                    rsq_data[all_th == False] = 0
                    alpha = (rsq_data - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
                    alpha[alpha > 1] = 1
                    alpha[alpha < 0] = 0

                    # Create flatmaps
                    print('Creating flatmaps...')
                    maps_names = []

                    # OD index (right - left: left eye stronger = red, right eye stronger = blue)
                    od_data = np.copy(od_mat[od_re_le_idx,...])
                    od_data[all_th == False] = np.nan
                    param_od = {'data': od_data,
                                'cmap': cmap_od,
                                'alpha': alpha,
                                'vmin': od_scale[0],
                                'vmax': od_scale[1],
                                'cbar': 'discrete',
                                'cortex_type': 'VertexRGB',
                                'description': 'Ocular dominance index (right - left)',
                                'curv_brightness': 0.6,
                                'curv_contrast': 0.2,
                                'add_roi': save_svg,
                                'cbar_label': 'OD (left eye <- -> right eye)',
                                'overlay_fn': overlay_fn,
                                'with_labels': True}
                    maps_names.append('od')

                    # OD z (right - left)
                    od_z_data = np.copy(od_mat[od_z_re_le_idx,...])
                    od_z_data[all_th == False] = np.nan
                    param_od_z = {'data': od_z_data,
                                  'cmap': cmap_od,
                                  'alpha': alpha,
                                  'vmin': od_z_scale[0],
                                  'vmax': od_z_scale[1],
                                  'cbar': 'discrete',
                                  'cortex_type': 'VertexRGB',
                                  'description': 'Ocular dominance z (right - left)',
                                  'curv_brightness': 0.6,
                                  'curv_contrast': 0.2,
                                  'add_roi': False,
                                  'cbar_label': 'OD z (left eye <- -> right eye)',
                                  'overlay_fn': overlay_fn,
                                  'with_labels': True}
                    maps_names.append('od_z')

                    # draw flatmaps
                    volumes = {}
                    for maps_name in maps_names:

                        # create flatmap
                        roi_name = 'prf_{}'.format(maps_name)
                        roi_param = {'subject': pycortex_subject,
                                     'xfmname': None,
                                     'roi_name': roi_name}
                        print(roi_name)
                        exec('param_{}.update(roi_param)'.format(maps_name))
                        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name=maps_name))
                        flatmap_fn = f'{flatmaps_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}-{maps_name}.pdf'
                        exec("plt.savefig('{}')".format(flatmap_fn))
                        print(f"Saved flatmap: {flatmap_fn}")
                        plt.close()

                        # save flatmap as dataset
                        exec('vol_description = param_{}["description"]'.format(maps_name))
                        exec('volume = volume_{}'.format(maps_name))
                        volumes.update({vol_description: volume})

                    # save dataset
                    dataset_file = f"{datasets_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}_od.hdf"
                    if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
                    dataset = cortex.Dataset(data=volumes)
                    dataset.save(dataset_file)
                    print(f"Saved dataset: {dataset_file}")
