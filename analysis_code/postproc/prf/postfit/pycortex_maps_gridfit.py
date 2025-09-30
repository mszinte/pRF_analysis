"""
-----------------------------------------------------------------------------------------
pycortex_maps_gridfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save in svg (e.g. no)
sys.argv[5]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_gridfit.py [main directory] [project name] 
                                   [subject num] [save_in_svg] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit/

python pycortex_maps_gridfit.py ~/disks/meso_S/data MotConf sub-01 n
python pycortex_maps_gridfit.py ~/disks/meso_S/data MotConf sub-170k n

python pycortex_maps_gridfit.py ~/disks/meso_S/data RetinoMaps sub-01 n
python pycortex_maps_gridfit.py ~/disks/meso_S/data RetinoMaps sub-170k n

python pycortex_maps_gridfit.py ~/disks/meso_S/data amblyo_prf sub-01 n
python pycortex_maps_gridfit.py ~/disks/meso_S/data amblyo_prf sub-170k n

python pycortex_maps_gridfit.py ~/disks/meso_S/data amblyo_prf sub-01 n prf_em_ctrl
python pycortex_maps_gridfit.py ~/disks/meso_S/data amblyo_prf sub-170k n prf_em_ctrl
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
import json
import cortex
import importlib
import numpy as np
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg_in = sys.argv[4]
if len(sys.argv) > 5: output_folder = sys.argv[5]
else: output_folder = "prf"
    
try:
    if save_svg_in == 'yes' or save_svg_in == 'y':
        save_svg = True
    elif save_svg_in == 'no' or save_svg_in == 'n':
        save_svg = False
    else:
        raise ValueError
except ValueError:
    sys.exit('Error: incorrect input (Yes, yes, y or No, no, n)')
if subject == 'sub-170k': save_svg = False
else: save_svg = save_svg

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']
maps_names_gauss = analysis_info['maps_names_gauss']

# Maps settings
for idx, col_name in enumerate(maps_names_gauss):
    exec("{}_idx = idx".format(col_name)) 
cmap_polar, cmap_uni, cmap_ecc_size = 'hsv', 'Reds', 'Spectral'
col_offset = 1.0/14.0
cmap_steps = 255

# plot scales
rsq_scale = [0, 1]
ecc_scale = [0, 7.5]
size_scale = [0, 7.5]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)
 
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    
    # define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/{}".format(main_dir, project_dir, 
                                                           subject, format_, output_folder)
    fit_dir = "{}/fit".format(prf_dir)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    flatmaps_dir = '{}/pycortex/flatmaps_avg_gauss_gridfit'.format(prf_dir)
    datasets_dir = '{}/pycortex/datasets_avg_gauss_gridfit'.format(prf_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    
    if format_ == 'fsnative':
        deriv_avg_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_prf-deriv_gauss_gridfit.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_avg_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_prf-deriv_gauss_gridfit.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        results = load_surface_pycortex(L_fn=deriv_avg_fn_L, 
                                        R_fn=deriv_avg_fn_R)
        deriv_mat = results['data_concat']
        
    elif format_ == '170k':
        deriv_avg_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-deriv_gauss_gridfit.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        results = load_surface_pycortex(brain_fn=deriv_avg_fn)
        deriv_mat = results['data_concat']
        if subject == 'sub-170k': save_svg = save_svg
        else: save_svg = False
    
    print('Creating flatmaps...')
    
    maps_names = []
    
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
    
    # r-square
    rsq_data = deriv_mat[prf_rsq_idx,...]
    alpha_range = analysis_info["alpha_range"]
    alpha = (rsq_data - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
    alpha[alpha>1]=1
    param_rsq = {'data': rsq_data, 
                 'cmap': cmap_uni, 
                 'alpha': alpha,
                 'vmin': rsq_scale[0], 
                 'vmax': rsq_scale[1], 
                 'cbar': 'discrete',
                 'cortex_type': 'VertexRGB',
                 'description': 'Gauss pRF R2',
                 'curv_brightness': 1,
                 'curv_contrast': 0.1, 
                 'add_roi': save_svg,
                 'cbar_label': 'pRF R2', 
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
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
        exec("plt.savefig('{}/{}_task-{}_{}.pdf')".format(
            flatmaps_dir, subject, prf_task_name, maps_name))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})
    
    # save dataset
    dataset_file = "{}/{}_task-{}_avg_gauss_gridfit.hdf".format(datasets_dir, subject, prf_task_name)
    if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)