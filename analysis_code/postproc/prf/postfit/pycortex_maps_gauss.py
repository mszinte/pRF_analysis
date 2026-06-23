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
sys.argv[4]: analysis name (e.g. prf)
sys.argv[5]: save in svg (e.g. no)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON PERSONAL SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_gauss.py [main directory] [project name] [analysis name]
                                 [subject num] [save_in_svg]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python pycortex_maps_gauss.py ~/disks/meso_S/data RetinoMaps sub-01 prf n
python pycortex_maps_gauss.py ~/disks/meso_S/data RetinoMaps sub-hcp1.6mm prf n
python pycortex_maps_gauss.py ~/disks/meso_S/data amblyo7T_prf sub-04 prf n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import glob
import cortex
import importlib
import numpy as np
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex


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
task_names = analysis_info['analysis_task_names']
maps_names_gauss = analysis_info['maps_names_gauss']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]

# plot scales
rsq_scale = analysis_info['flatmap_rsq_scale']
ecc_scale = analysis_info['flatmap_ecc_scale']
size_scale = analysis_info['flatmap_size_scale']
alpha_range = analysis_info["flatmap_alpha_range"]


# Maps settings
for idx, col_name in enumerate(maps_names_gauss):
    exec("{}_idx = idx".format(col_name)) 
cmap_polar, cmap_uni, cmap_ecc_size = 'hsv', 'Reds', 'Spectral'
col_offset = 1.0/14.0
cmap_steps = 255

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

for avg_method in avg_methods:
    if "loo" in avg_method:
        continue  # Skip if it contains "loo"

    
    for format_ in formats:

        rois_methods_format = rois_methods[format_]    
        for rois_method_format in rois_methods_format:
            if 'drawn' in rois_method_format:
                # no yet drawn rois as it is based on these results
                overlay_fn = 'overlays.svg'
                rois_method_format_txt = ""
            else: 
                overlay_fn = f"overlays_{rois_method_format}.svg"
                rois_method_format_txt = f"_{rois_method_format}"
        
            # define directories and fn
            prf_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{output_folder}"
            if not os.path.isdir(prf_dir):
                print(f"[SKIP] corr_dir not found for format={format_}: {prf_dir}")
                continue
            
            fit_dir = "{}/fit".format(prf_dir)
            prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
            flatmaps_dir = '{}/pycortex/flatmaps_gauss'.format(prf_dir)
            datasets_dir = '{}/pycortex/datasets_gauss'.format(prf_dir)
            os.makedirs(flatmaps_dir, exist_ok=True)
            os.makedirs(datasets_dir, exist_ok=True)
        
            for task_name in task_names:
                # get run numbers
                if "single-run" in avg_method:
                    if format_ == 'fsnative':
                        deriv_list = glob.glob(f'{prf_deriv_dir}/{subject}_task-{task_name}*_hemi-L_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-gauss{dm_name}_deriv.func.gii')
                    elif format_ == '170k':
                        deriv_list = glob.glob(f'{prf_deriv_dir}/{subject}_task-{task_name}*_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-gauss{dm_name}_deriv.dtseries.nii')
                    runs = sorted(set('_run-' + f.split('run-')[1].split('_')[0] for f in deriv_list))
                else:
                    runs = ['']
                    
                for run in runs:
                    if format_ == 'fsnative':
                        pycortex_subject = subject
                        deriv_avg_fn_L = f'{prf_deriv_dir}/{subject}_task-{task_name}{run}_hemi-L_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-gauss{dm_name}_deriv.func.gii'
                        deriv_avg_fn_R = f'{prf_deriv_dir}/{subject}_task-{task_name}{run}_hemi-R_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-gauss{dm_name}_deriv.func.gii'
                        results = load_surface_pycortex(L_fn=deriv_avg_fn_L, 
                                                        R_fn=deriv_avg_fn_R)
                        deriv_mat = results['data_concat']
                        
                    elif format_ == '170k':
                        pycortex_subject = pycortex_subject_template
                        deriv_avg_fn = f'{prf_deriv_dir}/{subject}_task-{task_name}{run}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-gauss{dm_name}_deriv.dtseries.nii'
                        results = load_surface_pycortex(brain_fn=deriv_avg_fn)
                        deriv_mat = results['data_concat']
                
                    # threshold data
                    deriv_mat_th = deriv_mat
                    amp_down =  deriv_mat_th[amplitude_idx,...] > 0
                    rsqr_th_down = deriv_mat_th[prf_rsq_idx,...] >= analysis_info['rsqr_th']
                    size_th_down = deriv_mat_th[prf_size_idx,...] >= analysis_info['size_th'][0]
                    size_th_up = deriv_mat_th[prf_size_idx,...] <= analysis_info['size_th'][1]
                    ecc_th_down = deriv_mat_th[prf_ecc_idx,...] >= analysis_info['ecc_th'][0]
                    ecc_th_up = deriv_mat_th[prf_ecc_idx,...] <= analysis_info['ecc_th'][1]
                    all_th = np.array((amp_down, rsqr_th_down, size_th_down, size_th_up, ecc_th_down, ecc_th_up)) 
                    deriv_mat[prf_rsq_idx,np.logical_and.reduce(all_th)==False]=0
        
                    # Create flatmaps
                    print('Creating flatmaps...')
                    maps_names = []
        
                    # r-square
                    rsq_data = deriv_mat[prf_rsq_idx,...]
                    alpha = (rsq_data - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
                    alpha[alpha>1]=1
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
                    pol_comp_num = deriv_mat[polar_real_idx,...] + 1j * deriv_mat[polar_imag_idx,...]
                    polar_ang = np.angle(pol_comp_num)
                    ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)
                    ang_norm = np.fmod(ang_norm + col_offset,1)
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
                    ecc_data = deriv_mat[prf_ecc_idx,...]
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
                    size_data = deriv_mat[prf_size_idx,...]
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
                    
                    # draw flatmaps
                    volumes = {}
                    for maps_name in maps_names:
                    
                        # create flatmap
                        roi_name = 'prf_{}'.format(maps_name)
                        roi_param = {'subject': pycortex_subject, 
                                     'roi_name': roi_name}
                        print(roi_name)
                        exec('param_{}.update(roi_param)'.format(maps_name))
                        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name=maps_name))
                        exec("plt.savefig('{}/{}_task-{}{}_{}_{}_{}_{}{}_gauss{}-{}.pdf')".format(flatmaps_dir, subject, analysis_name, run, preproc_prep, filtering, normalization, avg_method, rois_method_format_txt, dm_name, maps_name))
                        plt.close()
                    
                        # save flatmap as dataset
                        exec('vol_description = param_{}["description"]'.format(maps_name))
                        exec('volume = volume_{}'.format(maps_name))
                        volumes.update({vol_description:volume})
                    
                    # save dataset
                    dataset_file = f"{datasets_dir}/{subject}_task-{task_name}{run}_{preproc_prep}_{filtering}_{normalization}_{avg_method}{rois_method_format_txt}_{analysis_name}-gauss{dm_name}.hdf"
                    if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
                    dataset = cortex.Dataset(data=volumes)
                    dataset.save(dataset_file)