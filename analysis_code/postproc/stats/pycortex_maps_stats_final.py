"""
-----------------------------------------------------------------------------------------
pycortex_maps_stats_final.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape) 
1. cd to function
>> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/stats
2. run python command
>> python pycortex_maps_glm_final.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/postproc/stats
python pycortex_maps_stats_final.py ~/disks/meso_shared RetinoMaps sub-01 n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

#  Debug import 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json
import cortex
import numpy as np
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex, create_colormap

#Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']
tasks = analysis_info['task_final_stats']
alpha_range = analysis_info["alpha_range"]

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg_in = sys.argv[4]

try:
    if save_svg_in == 'yes' or save_svg_in == 'y':
        save_svg = True
    elif save_svg_in == 'no' or save_svg_in == 'n':
        save_svg = False
    else:
        raise ValueError
except ValueError:
    sys.exit('Error: incorrect input (Yes, yes, y or No, no, n)')
       
# Maps settings
all_idx, pur_idx, sac_idx, pur_sac_idx, prf_idx, prf_pur_idx, prf_sac_idx, \
        prf_pur_sac_idx = 0,1,2,3,4,5,6,7

rsq_idx = 2

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define/create colormap
colormap_name = 'stats_colors'
colormap_dict = {'n/a': (255, 255, 255),
                 'Pursuit':  (227, 119, 194),
                 'Saccade': (140, 86, 75),
                 'Pursuit_and_Saccade': (148, 103, 189),
                 'Vision': (214, 39, 40),
                 'Vision_and_Pursuit': (44, 160, 44), 
                 'Vision_and_Saccade': (255, 127, 14),
                 'Vision_and_Pursuit_and_Saccade': (31, 119, 180)}
create_colormap(cortex_dir=cortex_dir, 
                colormap_name=colormap_name, 
                colormap_dict=colormap_dict)


if subject == 'sub-170k':
    formats = ['170k'] 
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # Define directories and fn
    stats_dir = "{}/{}/derivatives/pp_data/{}/{}/final_stats".format(
        main_dir, project_dir, subject,format_)
    prf_stats_dir = "{}/{}/derivatives/pp_data/{}/{}/prf/prf_derivatives".format(
        main_dir, project_dir, subject,format_)
    glm_stats_dir = "{}/{}/derivatives/pp_data/{}/{}/glm/glm_derivatives".format(
        main_dir, project_dir,  subject,format_)
    stats_results_dir = "{}/results".format(stats_dir)
    
    flatmaps_dir = '{}/pycortex/flatmaps_stats'.format(stats_dir)
    datasets_dir = '{}/pycortex/datasets_stats'.format(stats_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    
    if format_ == 'fsnative': 
        deriv_stats_fn_L = '{}/{}_hemi-L_final-stats.func.gii'.format(
            stats_results_dir, subject)
        deriv_stats_fn_R = '{}/{}_hemi-R_final-stats.func.gii'.format(
            stats_results_dir, subject)
        
        prf_stats_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
            prf_stats_dir, subject, prf_task_name)
        prf_stats_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
            prf_stats_dir, subject, prf_task_name)
        sac_stats_fn_L = '{}/{}_task-SacLoc_hemi-L_fmriprep_dct_loo-avg_glm-stats.func.gii'.format(
            glm_stats_dir, subject)
        sac_stats_fn_R = '{}/{}_task-SacLoc_hemi-R_fmriprep_dct_loo-avg_glm-stats.func.gii'.format(
            glm_stats_dir, subject)
        pur_stats_fn_L = '{}/{}_task-PurLoc_hemi-L_fmriprep_dct_loo-avg_glm-stats.func.gii'.format(
            glm_stats_dir, subject)
        pur_stats_fn_R = '{}/{}_task-PurLoc_hemi-R_fmriprep_dct_loo-avg_glm-stats.func.gii'.format(
            glm_stats_dir, subject)
        
        #  Load data
        results_stats = load_surface_pycortex(L_fn=deriv_stats_fn_L, R_fn=deriv_stats_fn_R)
        final_mat = results_stats['data_concat']

        results_prf = load_surface_pycortex(L_fn=prf_stats_fn_L, R_fn=prf_stats_fn_R)
        prf_mat = results_prf['data_concat']

        results_sac = load_surface_pycortex(L_fn=sac_stats_fn_L, R_fn=sac_stats_fn_R)
        sac_mat = results_sac['data_concat']

        results_pur = load_surface_pycortex(L_fn=pur_stats_fn_L, R_fn=pur_stats_fn_R)
        pur_mat = results_pur['data_concat']
        
    elif format_ == '170k':
        stats_avg_fn = '{}/{}_final-stats.dtseries.nii'.format(
            stats_results_dir, subject)
        prf_stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_loo-avg_prf-stats.dtseries.nii'.format(
            prf_stats_dir, subject, prf_task_name)
        sac_stats_avg_fn = '{}/{}_task-SacLoc_fmriprep_dct_loo-avg_glm-stats.dtseries.nii'.format(
            glm_stats_dir, subject)
        pur_stats_avg_fn = '{}/{}_task-PurLoc_fmriprep_dct_loo-avg_glm-stats.dtseries.nii'.format(
            glm_stats_dir, subject)
        
        #  Load data
        results_stats = load_surface_pycortex(brain_fn=stats_avg_fn)
        final_mat = results_stats['data_concat']

        results_prf = load_surface_pycortex(brain_fn=prf_stats_avg_fn)
        prf_mat = results_prf['data_concat']
        
        results_sac = load_surface_pycortex(brain_fn=sac_stats_avg_fn)
        sac_mat = results_sac['data_concat']

        results_pur = load_surface_pycortex(brain_fn=pur_stats_avg_fn)
        pur_mat = results_pur['data_concat']

        if subject == 'sub-170k':
            save_svg = save_svg
        else: 
            save_svg = False
    
    # Compute R2 from R
    prf_mat[rsq_idx,:] = prf_mat[rsq_idx,:]**2
    sac_mat[rsq_idx,:] = sac_mat[rsq_idx,:]**2
    pur_mat[rsq_idx,:] = pur_mat[rsq_idx,:]**2

    # threshold data
    # pRF
    prf_mat_th = prf_mat
    rsqr_th_down = prf_mat_th[rsq_idx,...] >= analysis_info['rsqr_th'][0]
    rsqr_th_up = prf_mat_th[rsq_idx,...] <= analysis_info['rsqr_th'][1]

    all_th = np.array((rsqr_th_down, rsqr_th_up)) 
    prf_mat[rsq_idx,np.logical_and.reduce(all_th)==False]=0
    
    # SacLoc
    sac_mat_th = sac_mat
    rsqr_th_down = sac_mat_th[rsq_idx,...] >= analysis_info['rsqr_th'][0]
    rsqr_th_up = sac_mat_th[rsq_idx,...] <= analysis_info['rsqr_th'][1]

    all_th = np.array((rsqr_th_down, rsqr_th_up)) 
    sac_mat[rsq_idx,np.logical_and.reduce(all_th)==False]=0
    
    # PurLoc
    pur_mat_th = pur_mat
    rsqr_th_down = pur_mat_th[rsq_idx,...] >= analysis_info['rsqr_th'][0]
    rsqr_th_up = pur_mat_th[rsq_idx,...] <= analysis_info['rsqr_th'][1]

    all_th = np.array((rsqr_th_down, rsqr_th_up)) 
    pur_mat[rsq_idx,np.logical_and.reduce(all_th)==False]=0


    #  Creat R2 for the all flatmap
    rsq_all = np.zeros((final_mat[all_idx,...].shape))
    for vert, categorie in enumerate(final_mat[all_idx,...]):
        if categorie == 1: rsq_all[vert] = pur_mat[rsq_idx,vert]
        
        elif categorie == 2: rsq_all[vert] = sac_mat[rsq_idx,vert]
        
        elif categorie == 3: rsq_all[vert] = np.nanmean([pur_mat[rsq_idx,vert], 
                                                         sac_mat[rsq_idx,vert]],
                                                        axis=0)
        
        elif categorie == 4: rsq_all[vert] = prf_mat[rsq_idx,vert]
        
        elif categorie == 5: rsq_all[vert] = np.nanmean([prf_mat[rsq_idx,vert], 
                                                         pur_mat[rsq_idx,vert]],
                                                        axis=0)
        
        elif categorie == 6: rsq_all[vert] = np.nanmean([prf_mat[rsq_idx,vert], 
                                                         sac_mat[rsq_idx,vert]],
                                                        axis=0)
        
        elif categorie == 7: rsq_all[vert] = np.nanmean([prf_mat[rsq_idx,vert], 
                                                         sac_mat[rsq_idx,vert], 
                                                         pur_mat[rsq_idx,vert]], 
                                                        axis=0)
            
    #  Creat R2 for the pur_sac flatmap       
    rsq_pur_sac = np.zeros((final_mat[pur_sac_idx,...].shape))
    for vert, categorie in enumerate(final_mat[pur_sac_idx,...]):
        if categorie == 3: rsq_pur_sac[vert] = np.nanmean([pur_mat[rsq_idx,vert], 
                                                           sac_mat[rsq_idx,vert]], 
                                                          axis=0)
            
    #  Creat R2 for the prf_pur flatmap               
    rsq_prf_pur = np.zeros((final_mat[prf_pur_idx,...].shape))
    for vert, categorie in enumerate(final_mat[prf_pur_idx,...]):
        if categorie == 5: rsq_prf_pur[vert] = np.nanmean([pur_mat[rsq_idx,vert], 
                                                           sac_mat[rsq_idx,vert]], 
                                                          axis=0)
            
    #  Creat R2 for the prf_sac flatmap               
    rsq_prf_sac = np.zeros((final_mat[prf_sac_idx,...].shape))
    for vert, categorie in enumerate(final_mat[prf_sac_idx,...]):
        if categorie == 6: rsq_prf_sac[vert] = np.nanmean([prf_mat[rsq_idx,vert], 
                                                           sac_mat[rsq_idx,vert]], 
                                                          axis=0)
    
    #  Creat R2 for the prf_pur_sac flatmap               
    rsq_prf_pur_sac = np.zeros((final_mat[all_idx,...].shape))
    for vert, categorie in enumerate(final_mat[all_idx,...]):
        if categorie == 7: rsq_prf_pur_sac[vert] = np.nanmean([prf_mat[rsq_idx,vert], 
                                                               sac_mat[rsq_idx,vert], 
                                                               pur_mat[rsq_idx,vert]], 
                                                              axis=0)
        
    # Create flatmaps  
    print('Creating flatmaps...')
    maps_names = []        
    
    #  Creat the all flatmap
    alpha_all = (rsq_all - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_all[alpha_all>1]=1
    
    final_data_all = final_mat[all_idx,...]
    param_all = {'data': final_data_all, 
                 'cmap': colormap_name, 
                 'alpha': alpha_all, 
                 'vmin': 0, 
                 'vmax': 7, 
                 'cbar': 'discrete_personalized', 
                 'cmap_steps': len(colormap_dict),
                 'cmap_dict': colormap_dict,
                 'cortex_type': 'VertexRGB', 
                 'description': 'final map', 
                 'curv_brightness': 0.1, 
                 'curv_contrast': 0.25, 
                 'add_roi': save_svg, 
                 'cbar_label': '', 
                 'with_labels': True}
    maps_names.append('all')
    
    #  Creat the pursuit flatmap
    rsq_pur = pur_mat[rsq_idx, :]
    alpha_pur = (rsq_pur - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_pur[alpha_pur>1]=1

    final_data_pur = final_mat[pur_idx,...]
    param_pursuit = {'data': final_data_pur, 
                     'cmap': colormap_name, 'alpha': alpha_pur, 
                     'vmin': 0, 
                     'vmax': 7, 
                     'cbar': 'discrete_personalized', 
                     'cmap_steps': len(colormap_dict),
                     'cmap_dict': colormap_dict,
                     'cortex_type': 'VertexRGB', 
                     'description': 'final map',
                     'curv_brightness': 0.1, 
                     'curv_contrast': 0.25, 
                     'add_roi': save_svg,
                     'cbar_label': '', 
                     'with_labels': True}
    maps_names.append('pursuit')
    
    #  Creat the saccade flatmap
    rsq_sac = sac_mat[rsq_idx, :]
    alpha_sac = (rsq_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_sac[alpha_sac>1]=1
    
    final_data_sac = final_mat[sac_idx,...]
    param_saccade = {'data': final_data_sac, 
                     'cmap': colormap_name, 'alpha': alpha_sac, 
                     'vmin': 0, 
                     'vmax': 7, 
                     'cbar': 'discrete_personalized', 
                     'cmap_steps': len(colormap_dict),
                     'cmap_dict': colormap_dict,
                     'cortex_type': 'VertexRGB', 
                     'description': 'final map',
                     'curv_brightness': 0.1, 
                     'curv_contrast': 0.25, 
                     'add_roi': save_svg,
                     'cbar_label': '', 
                     'with_labels': True}
    maps_names.append('saccade')
    
    #  Creat the pursuit_and_saccade flatmap
    alpha_pur_sac = (rsq_pur_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_pur_sac[alpha_pur_sac>1]=1
    
    final_data_pur_sac = final_mat[pur_sac_idx,...]
    param_pursuit_and_saccade = {'data': final_data_pur_sac, 
                                 'cmap': colormap_name, 
                                 'alpha': alpha_pur_sac, 
                                 'vmin': 0, 
                                 'vmax': 7, 
                                 'cbar': 'discrete_personalized', 
                                 'cmap_steps': len(colormap_dict),
                                 'cmap_dict': colormap_dict,
                                 'cortex_type': 'VertexRGB', 
                                 'description': 'final map',
                                 'curv_brightness': 0.1, 
                                 'curv_contrast': 0.25, 
                                 'add_roi': save_svg,
                                 'cbar_label': '', 
                                 'with_labels': True}
    maps_names.append('pursuit_and_saccade')
    
    #  Creat the vision flatmap
    rsq_prf = prf_mat[rsq_idx, :]
    alpha_prf = (rsq_prf - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_prf[alpha_prf>1]=1

    final_data_prf = final_mat[prf_idx,...]
    param_vision = {'data': final_data_prf, 
                    'cmap': colormap_name, 
                    'alpha': alpha_prf, 
                    'vmin': 0, 
                    'vmax': 7, 
                    'cbar': 'discrete_personalized', 
                    'cmap_steps': len(colormap_dict),
                    'cmap_dict': colormap_dict,
                    'cortex_type': 'VertexRGB', 
                    'description': 'final map',
                    'curv_brightness': 0.1, 
                    'curv_contrast': 0.25, 
                    'add_roi': save_svg,
                    'cbar_label': '',
                    'with_labels': True}
    maps_names.append('vision')
    
    #  Creat the vision_and_pursuite flatmap
    alpha_prf_pur = (rsq_prf_pur - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_prf_pur[alpha_prf_pur>1]=1
    
    final_data_prf_pur = final_mat[prf_pur_idx,...]
    param_vision_and_pursuit = {'data': final_data_prf_pur, 
                                'cmap': colormap_name, 
                                'alpha': alpha_prf_pur, 
                                'vmin': 0, 
                                'vmax': 7, 
                                'cbar': 'discrete_personalized', 
                                'cmap_steps': len(colormap_dict), 
                                'cmap_dict': colormap_dict,
                                'cortex_type': 'VertexRGB', 
                                'description': 'final map',
                                'curv_brightness': 0.1, 
                                'curv_contrast': 0.25, 
                                'add_roi': save_svg,
                                'cbar_label': '', 
                                'with_labels': True}
    maps_names.append('vision_and_pursuit')
    
    #  Creat the vision_and_saccade flatmap
    alpha_prf_sac = (rsq_prf_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_prf_sac[alpha_prf_sac>1]=1
    
    final_data_prf_sac = final_mat[prf_sac_idx,...]
    param_vision_and_saccade = {'data': final_data_prf_sac, 
                                'cmap': colormap_name, 
                                'alpha': alpha_prf_sac, 
                                'vmin': 0, 
                                'vmax': 7, 
                                'cbar': 'discrete_personalized', 
                                'cmap_steps': len(colormap_dict), 
                                'cmap_dict': colormap_dict,
                                'cortex_type': 'VertexRGB', 
                                'description': 'final map',
                                'curv_brightness': 0.1, 
                                'curv_contrast': 0.25, 
                                'add_roi': save_svg,
                                'cbar_label': '', 
                                'with_labels': True}
    maps_names.append('vision_and_saccade')
    
    #  Creat the vision_and_pursuit_and_saccade flatmap
    alpha_prf_pur_sac = (rsq_prf_pur_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha_prf_pur_sac[alpha_prf_pur_sac>1]=1

    final_data_prf_pur_sac = final_mat[prf_pur_sac_idx,...]
    param_vision_and_pursuit_and_saccade = {'data': final_data_prf_pur_sac, 
                                            'cmap': colormap_name, 
                                            'alpha': alpha_prf_pur_sac, 
                                            'vmin': 0, 'vmax': 7, 
                                            'cbar': 'discrete_personalized', 
                                            'cmap_steps': len(colormap_dict), 
                                            'cmap_dict': colormap_dict,
                                            'cortex_type': 'VertexRGB', 
                                            'description': 'final map', 
                                            'curv_brightness': 0.1, 
                                            'curv_contrast': 0.25, 
                                            'add_roi': save_svg, 
                                            'cbar_label': '',  
                                            'with_labels': True}
    maps_names.append('vision_and_pursuit_and_saccade')

    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = 'final_stats_{}'.format(maps_name)
        roi_param = {'subject': pycortex_subject, 'xfmname': None, 'roi_name': roi_name}
        print(roi_name)
        exec('param_{}.update(roi_param)'.format(maps_name))
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
        exec("plt.savefig('{}/{}_final_{}.pdf')".format(flatmaps_dir, subject, maps_name))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})
    
    # save dataset
    dataset_file = "{}/{}_final.hdf".format(datasets_dir, subject)
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)
    
    
