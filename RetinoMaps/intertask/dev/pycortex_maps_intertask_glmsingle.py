"""
-----------------------------------------------------------------------------------------
pycortex_maps_intertask_glmsingle.py.py
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
>> cd ~/disks/meso_H/projects/pRF_analysis/RetinoMaps/intertask/dev
2. run python command
>> python pycortex_maps_intertask.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/RetinoMaps/intertask/dev
python pycortex_maps_intertask_glmsingle.py ~/disks/meso_shared RetinoMaps sub-01 n
python pycortex_maps_intertask_glmsingle.py ~/disks/meso_shared RetinoMaps sub-170k n
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
import copy
import json
import cortex
import numpy as np
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex, create_colormap

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg_in = sys.argv[4]

#Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']
alpha_range = analysis_info["alpha_range"]
group_tasks = analysis_info['task_intertask']

# group_tasks = [['PurLoc', 'SacLoc', 'pRF']] # To remove !!!!

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

slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define/create colormap
# colormap_name = 'stats_colors'
colormap_name = 'viridis'
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

for tasks in group_tasks : 
    if 'SacVELoc' in tasks: 
        print(tasks)
        suffix = 'SacVE_PurVE'
        sac_task = 'SacVELoc'
        pur_task = 'PurVELoc'
        print(suffix)
    else : 
        print(tasks)
        suffix = 'Sac_Pur'
        sac_task = 'SacLoc'
        pur_task = 'PurLoc'
        print(suffix)

    for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
        # Define directories and fn
        intertask_dir = "{}/{}/derivatives/pp_data/{}/{}/intertask_glmsingle".format(
            main_dir, project_dir, subject,format_)
        prf_stats_dir = "{}/{}/derivatives/pp_data/{}/{}/prf/prf_derivatives".format(
            main_dir, project_dir, subject,format_)
        glm_stats_dir = "{}/{}/derivatives/pp_data/{}/{}/glm/glm_single_derivatives".format(
            main_dir, project_dir,  subject,format_)
        
        flatmaps_dir = '{}/pycortex/flatmaps_stats'.format(intertask_dir)
        datasets_dir = '{}/pycortex/datasets_stats'.format(intertask_dir)
        
        os.makedirs(flatmaps_dir, exist_ok=True)
        os.makedirs(datasets_dir, exist_ok=True)
        
        if format_ == 'fsnative': 
            deriv_stats_fn_L = '{}/{}_hemi-L_intertask_{}.func.gii'.format(
                intertask_dir, subject, suffix)
            deriv_stats_fn_R = '{}/{}_hemi-R_intertask_{}.func.gii'.format(
                intertask_dir, subject, suffix)
            
            prf_stats_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_prf-stats_loo-median.func.gii'.format(
                prf_stats_dir, subject, prf_task_name)
            prf_stats_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_prf-stats_loo-median.func.gii'.format(
                prf_stats_dir, subject, prf_task_name)
            
            sac_stats_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_glm-stats_loo-median.func.gii'.format(
                glm_stats_dir, subject, sac_task)
            sac_stats_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_glm-stats_loo-median.func.gii'.format(
                glm_stats_dir, subject, sac_task)
            
            pur_stats_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_glm-stats_loo-median.func.gii'.format(
                glm_stats_dir, subject, pur_task)
            pur_stats_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_glm-stats_loo-median.func.gii'.format(
                glm_stats_dir, subject, pur_task)
            
            #  Load data
            results_stats = load_surface_pycortex(L_fn=deriv_stats_fn_L, R_fn=deriv_stats_fn_R)
            intertask_mat = results_stats['data_concat']
    
            results_prf = load_surface_pycortex(L_fn=prf_stats_fn_L, R_fn=prf_stats_fn_R)
            prf_mat = results_prf['data_concat']
    
            results_sac = load_surface_pycortex(L_fn=sac_stats_fn_L, R_fn=sac_stats_fn_R)
            sac_mat = results_sac['data_concat']
    
            results_pur = load_surface_pycortex(L_fn=pur_stats_fn_L, R_fn=pur_stats_fn_R)
            pur_mat = results_pur['data_concat']
            
        elif format_ == '170k':
            stats_avg_fn = '{}/{}_intertask_{}.dtseries.nii'.format(
                intertask_dir, subject, suffix)
            prf_stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii'.format(
                prf_stats_dir, subject, prf_task_name)
            sac_stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_avg_glm-stats_loo-median.dtseries.nii'.format(
                glm_stats_dir, subject, sac_task)
            pur_stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_avg_glm-stats_loo-median.dtseries.nii'.format(
                glm_stats_dir, subject, pur_task)
            
            #  Load data
            results_stats = load_surface_pycortex(brain_fn=stats_avg_fn)
            intertask_mat = results_stats['data_concat']
    
            results_prf = load_surface_pycortex(brain_fn=prf_stats_avg_fn)
            prf_mat = results_prf['data_concat']
            
            results_sac = load_surface_pycortex(brain_fn=sac_stats_avg_fn)
            sac_mat = results_sac['data_concat']
    
            results_pur = load_surface_pycortex(brain_fn=pur_stats_avg_fn)
            pur_mat = results_pur['data_concat']
    
            if subject == 'sub-170k': save_svg = save_svg
            else: save_svg = False
        

        # Compute R2 from R
        prf_mat[rvalue_idx,:] =  prf_mat[rvalue_idx,:]**2
        sac_mat[rvalue_idx,:] = sac_mat[rvalue_idx,:]**2
        pur_mat[rvalue_idx,:] = pur_mat[rvalue_idx,:]**2
    
        # threshold data
        # pRF
        prf_mat_corrected = copy.copy(prf_mat)
        prf_mat_corrected_th = prf_mat_corrected
        if analysis_info['stats_th'] == 0.05: stats_th_down = prf_mat_corrected_th[corr_pvalue_5pt_idx,...] <= 0.05
        elif analysis_info['stats_th'] == 0.01: stats_th_down = prf_mat_corrected_th[corr_pvalue_1pt_idx,...] <= 0.01
        prf_mat_corrected[rvalue_idx, stats_th_down==False]=0 # put this to zero to not plot it

        # SacLoc        
        sac_mat_corrected = copy.copy(sac_mat)
        sac_mat_corrected_th = sac_mat_corrected
        if analysis_info['stats_th'] == 0.05: stats_th_down = sac_mat_corrected_th[corr_pvalue_5pt_idx,...] <= 0.05
        elif analysis_info['stats_th'] == 0.01: stats_th_down = sac_mat_corrected_th[corr_pvalue_1pt_idx,...] <= 0.01
        sac_mat_corrected[rvalue_idx, stats_th_down==False]=0 # put this to zero to not plot it
 
        # PurLoc        
        pur_mat_corrected = copy.copy(pur_mat)
        pur_mat_corrected_th = pur_mat_corrected
        if analysis_info['stats_th'] == 0.05: stats_th_down = pur_mat_corrected_th[corr_pvalue_5pt_idx,...] <= 0.05
        elif analysis_info['stats_th'] == 0.01: stats_th_down = pur_mat_corrected_th[corr_pvalue_1pt_idx,...] <= 0.01
        pur_mat_corrected[rvalue_idx, stats_th_down==False]=0 # put this to zero to not plot it
        
        
        #  Creat R2 for the all flatmap
        rsq_all = np.zeros((intertask_mat[all_idx,...].shape))
        for vert, categorie in enumerate(intertask_mat[all_idx,...]):
            if categorie == 1: rsq_all[vert] = pur_mat_corrected[rvalue_idx,vert]
            
            elif categorie == 2: rsq_all[vert] = sac_mat_corrected[rvalue_idx,vert]
            
            elif categorie == 3: rsq_all[vert] = np.nanmedian([pur_mat_corrected[rvalue_idx,vert], 
                                                             sac_mat_corrected[rvalue_idx,vert]],
                                                            axis=0)
            
            elif categorie == 4: rsq_all[vert] = prf_mat_corrected[rvalue_idx,vert]
            
            elif categorie == 5: rsq_all[vert] = np.nanmedian([prf_mat_corrected[rvalue_idx,vert], 
                                                             pur_mat_corrected[rvalue_idx,vert]],
                                                            axis=0)
            
            elif categorie == 6: rsq_all[vert] = np.nanmedian([prf_mat_corrected[rvalue_idx,vert], 
                                                             sac_mat_corrected[rvalue_idx,vert]],
                                                            axis=0)
            
            elif categorie == 7: rsq_all[vert] = np.nanmedian([prf_mat_corrected[rvalue_idx,vert], 
                                                              sac_mat_corrected[rvalue_idx,vert], 
                                                              pur_mat_corrected[rvalue_idx,vert]], 
                                                            axis=0)
                
        #  Creat R2 for the pur_sac flatmap       
        rsq_pur_sac = np.zeros((intertask_mat[pur_sac_idx,...].shape))
        for vert, categorie in enumerate(intertask_mat[pur_sac_idx,...]):
            if categorie == 3: rsq_pur_sac[vert] = np.nanmedian([pur_mat_corrected[rvalue_idx,vert], 
                                                               sac_mat_corrected[rvalue_idx,vert]], 
                                                              axis=0)
                
        #  Creat R2 for the prf_pur flatmap               
        rsq_prf_pur = np.zeros((intertask_mat[prf_pur_idx,...].shape))
        for vert, categorie in enumerate(intertask_mat[prf_pur_idx,...]):
            if categorie == 5: rsq_prf_pur[vert] = np.nanmedian([pur_mat_corrected[rvalue_idx,vert], 
                                                               sac_mat_corrected[rvalue_idx,vert]], 
                                                              axis=0)
                
        #  Creat R2 for the prf_sac flatmap               
        rsq_prf_sac = np.zeros((intertask_mat[prf_sac_idx,...].shape))
        for vert, categorie in enumerate(intertask_mat[prf_sac_idx,...]):
            if categorie == 6: rsq_prf_sac[vert] = np.nanmedian([prf_mat_corrected[rvalue_idx,vert], 
                                                               sac_mat_corrected[rvalue_idx,vert]], 
                                                              axis=0)
        
        #  Creat R2 for the prf_pur_sac flatmap               
        rsq_prf_pur_sac = np.zeros((intertask_mat[all_idx,...].shape))
        for vert, categorie in enumerate(intertask_mat[all_idx,...]):
            if categorie == 7: rsq_prf_pur_sac[vert] = np.nanmean([prf_mat[rvalue_idx,vert], 
                                                                   sac_mat_corrected[rvalue_idx,vert], 
                                                                   pur_mat_corrected[rvalue_idx,vert]], 
                                                                  axis=0)
            
        # Create flatmaps  
        print('Creating flatmaps...')
        maps_names = []        
        
        #  Creat the all flatmap
        alpha_all = (rsq_all - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_all[alpha_all>1]=1
        
        
        intertask_data_all = intertask_mat[all_idx,...]
        alpha_all[intertask_data_all == 0] = 0
        param_all = {'data': rsq_all, 
                     'cmap': colormap_name, 
                     'alpha': alpha_all, 
                     'vmin': 0, 
                     'vmax': 1, 
                     # 'cbar': 'discrete_personalized', 
                     'cmap_steps': len(colormap_dict),
                     'cmap_dict': colormap_dict,
                     'cortex_type': 'VertexRGB', 
                     'description': 'intertask map', 
                     'curv_brightness': 0.1, 
                     'curv_contrast': 0.25, 
                     'add_roi': save_svg, 
                     'cbar_label': '', 
                     'with_labels': True}
        maps_names.append('all')
        
        #  Creat the pursuit flatmap
        rsq_pur = pur_mat_corrected[rvalue_idx, :]
        
        
        alpha_pur = (rsq_pur - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_pur[alpha_pur>1]=1

        
    
        intertask_data_pur = intertask_mat[pur_idx,...]
        alpha_pur[intertask_data_pur == 0] = 0
        param_pursuit = {'data': rsq_pur, 
                         'cmap': colormap_name, 
                         'alpha': alpha_pur, 
                         'vmin': 0, 
                         'vmax': 1, 
                         # 'cbar': 'discrete_personalized', 
                         'cmap_steps': len(colormap_dict),
                         'cmap_dict': colormap_dict,
                         'cortex_type': 'VertexRGB', 
                         'description': 'intertask map',
                         'curv_brightness': 0.1, 
                         'curv_contrast': 0.25, 
                         'add_roi': save_svg,
                         'cbar_label': '', 
                         'with_labels': True}
        maps_names.append('pursuit')
        
        #  Creat the saccade flatmap
        rsq_sac = sac_mat_corrected[rvalue_idx, :]
        alpha_sac = (rsq_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_sac[alpha_sac>1]=1
        
        
        intertask_data_sac = intertask_mat[sac_idx,...]
        alpha_sac[intertask_data_sac == 0] = 0
        param_saccade = {'data': rsq_sac, 
                         'cmap': colormap_name, 
                         'alpha': alpha_sac, 
                         'vmin': 0, 
                         'vmax': 1, 
                         # 'cbar': 'discrete_personalized', 
                         'cmap_steps': len(colormap_dict),
                         'cmap_dict': colormap_dict,
                         'cortex_type': 'VertexRGB', 
                         'description': 'intertask map',
                         'curv_brightness': 0.1, 
                         'curv_contrast': 0.25, 
                         'add_roi': save_svg,
                         'cbar_label': '', 
                         'with_labels': True}
        maps_names.append('saccade')
        
        #  Creat the pursuit_and_saccade flatmap
        alpha_pur_sac = (rsq_pur_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_pur_sac[alpha_pur_sac>1]=1
        
        
        intertask_data_pur_sac = intertask_mat[pur_sac_idx,...]
        alpha_pur_sac[intertask_data_pur_sac == 0] = 0
        param_pursuit_and_saccade = {'data': rsq_pur_sac, 
                                     'cmap': colormap_name, 
                                     'alpha': alpha_pur_sac, 
                                     'vmin': 0, 
                                     'vmax': 1, 
                                     # 'cbar': 'discrete_personalized', 
                                     'cmap_steps': len(colormap_dict),
                                     'cmap_dict': colormap_dict,
                                     'cortex_type': 'VertexRGB', 
                                     'description': 'intertask map',
                                     'curv_brightness': 0.1, 
                                     'curv_contrast': 0.25, 
                                     'add_roi': save_svg,
                                     'cbar_label': '', 
                                     'with_labels': True}
        maps_names.append('pursuit_and_saccade')
        
        #  Creat the vision flatmap
        rsq_prf = prf_mat_corrected[rvalue_idx, :].astype(np.float64)
        alpha_prf = (rsq_prf - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_prf[alpha_prf>1]=1
        
    
        intertask_data_prf = intertask_mat[prf_idx,...]
        alpha_prf[intertask_data_prf == 0] = 0
        param_vision = {'data': rsq_prf, 
                        'cmap': colormap_name, 
                        'alpha': alpha_prf, 
                        'vmin': 0, 
                        'vmax': 1, 
                        # 'cbar': 'discrete_personalized', 
                        'cmap_steps': len(colormap_dict),
                        'cmap_dict': colormap_dict,
                        'cortex_type': 'VertexRGB', 
                        'description': 'intertask map',
                        'curv_brightness': 0.1, 
                        'curv_contrast': 0.25, 
                        'add_roi': save_svg,
                        'cbar_label': '',
                        'with_labels': True}
        maps_names.append('vision')
        
        #  Creat the vision_and_pursuite flatmap
        alpha_prf_pur = (rsq_prf_pur - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_prf_pur[alpha_prf_pur>1]=1
        
        intertask_data_prf_pur = intertask_mat[prf_pur_idx,...]
        alpha_prf_pur[intertask_data_prf_pur == 0] = 0
        param_vision_and_pursuit = {'data': rsq_prf_pur, 
                                    'cmap': colormap_name, 
                                    'alpha': alpha_prf_pur, 
                                    'vmin': 0, 
                                    'vmax': 1, 
                                    # 'cbar': 'discrete_personalized', 
                                    'cmap_steps': len(colormap_dict), 
                                    'cmap_dict': colormap_dict,
                                    'cortex_type': 'VertexRGB', 
                                    'description': 'intertask map',
                                    'curv_brightness': 0.1, 
                                    'curv_contrast': 0.25, 
                                    'add_roi': save_svg,
                                    'cbar_label': '', 
                                    'with_labels': True}
        maps_names.append('vision_and_pursuit')
        
        #  Creat the vision_and_saccade flatmap
        alpha_prf_sac = (rsq_prf_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_prf_sac[alpha_prf_sac>1]=1
        
        
        intertask_data_prf_sac = intertask_mat[prf_sac_idx,...]
        alpha_prf_sac[intertask_data_prf_sac == 0] = 0
        param_vision_and_saccade = {'data': rsq_prf_sac, 
                                    'cmap': colormap_name, 
                                    'alpha': alpha_prf_sac, 
                                    'vmin': 0, 
                                    'vmax': 1, 
                                    # 'cbar': 'discrete_personalized', 
                                    'cmap_steps': len(colormap_dict), 
                                    'cmap_dict': colormap_dict,
                                    'cortex_type': 'VertexRGB', 
                                    'description': 'intertask map',
                                    'curv_brightness': 0.1, 
                                    'curv_contrast': 0.25, 
                                    'add_roi': save_svg,
                                    'cbar_label': '', 
                                    'with_labels': True}
        maps_names.append('vision_and_saccade')
        
        #  Creat the vision_and_pursuit_and_saccade flatmap
        alpha_prf_pur_sac = (rsq_prf_pur_sac - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha_prf_pur_sac[alpha_prf_pur_sac>1]=1
        
    
        intertask_data_prf_pur_sac = intertask_mat[prf_pur_sac_idx,...]
        alpha_prf_pur_sac[intertask_data_prf_pur_sac == 0] = 0
        param_vision_and_pursuit_and_saccade = {'data': rsq_prf_pur_sac, 
                                                'cmap': colormap_name, 
                                                'alpha': alpha_prf_pur_sac, 
                                                'vmin': 0, 
                                                'vmax': 1, 
                                                # 'cbar': 'discrete_personalized', 
                                                'cmap_steps': len(colormap_dict), 
                                                'cmap_dict': colormap_dict,
                                                'cortex_type': 'VertexRGB', 
                                                'description': 'intertask map', 
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
            roi_name = '{}'.format(maps_name)
            roi_param = {'subject': pycortex_subject, 'xfmname': None, 'roi_name': roi_name}
            print(roi_name)
            exec('param_{}.update(roi_param)'.format(maps_name))
            exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
            exec("plt.savefig('{}/{}_{}_{}.pdf')".format(flatmaps_dir, subject, maps_name, suffix))
            plt.close()
        
            # save flatmap as dataset
            exec('vol_description = param_{}["description"]'.format(maps_name))
            exec('volume = volume_{}'.format(maps_name))
            volumes.update({vol_description:volume})
        
        # save dataset
        dataset_file = "{}/{}_{}_{}.hdf".format(datasets_dir, subject, maps_name, suffix)
        dataset = cortex.Dataset(data=volumes)
        dataset.save(dataset_file)
    
    
