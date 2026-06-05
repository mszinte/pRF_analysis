"""
-----------------------------------------------------------------------------------------
pycortex_maps_ncsf.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset of nCSF results
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save in svg (e.g. no)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
2. run python command
>> python pycortex_maps_ncsf.py [main directory] [project] [subject] [save_svg_in] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
python pycortex_maps_ncsf.py ~/disks/meso_shared nCSF sub-01 n
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net) 
and Martin Szinte (martin.szinte@gmail.com)
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
import numpy as np
import matplotlib.pyplot as plt


# Personal import
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex, create_colormap

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
nCSF_settings_path = os.path.join(base_dir, project_dir, "nCSF-analysis.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, nCSF_settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

# formats = analysis_info['formats']
# extensions = analysis_info['extensions']

formats = [analysis_info['formats'][0]] #######################################
extensions = [analysis_info['extensions'][0]]

ncsf_task_name = analysis_info['nCSF_task_name']
maps_names_ncsf = analysis_info['maps_names_ncsf']
maps_names_ncsf_stats = analysis_info['maps_names_ncsf_stats']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# plot scales
rsq_scale = analysis_info['flatmap_ncsf_rsq_scale']
SFp_scale = analysis_info['flatmap_ncsf_SFp_scale']
CSp_scale = analysis_info['flatmap_ncsf_CSp_scale']
auc_scale = analysis_info['flatmap_ncsf_auc_scale']
crf_exp_scale = analysis_info['flatmap_ncsf_crf_exp_scale']
normalize_auc_scale = analysis_info['flatmap_ncsf_normalize_auc_scale']
SFmax_scale = analysis_info['flatmap_ncsf_SFmax_scale']
alpha_range = analysis_info["flatmap_alpha_range"]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

cmap_SFp_CSp_SFmax, cmap_r2, cmap_auc = 'viridis', 'Reds', 'viridis'
col_offset = 1.0/14.0
cmap_steps = 255

for avg_method in avg_methods:
    # define maps name
    if 'loo' in avg_method: maps_names_ncsf = maps_names_ncsf + ['ncsf_loo_rsq']
    else: maps_names_css = maps_names_ncsf
    
    for idx, col_name in enumerate(maps_names_ncsf + maps_names_ncsf_stats ):
        exec("{}_idx = idx".format(col_name))
        
    # define rsquared or loo-rsquared idx
    if 'loo' in avg_method:
        rsq_idx2use = ncsf_loo_rsq_idx
        rsq_description = 'nCSF pRF LOO R2'
        rsq_cbar_label = 'nCSF LOO R2'
    else:
        rsq_idx2use = ncsf_rsq_idx
        rsq_description = 'CSS pRF R2'
        rsq_cbar_label = 'pRF R2'  
        
    for format_ in formats:

        rois_methods_format = rois_methods[format_]    
        for rois_method_format in rois_methods_format:
            overlay_fn = f"overlays_{rois_method_format}.svg"
            rois_method_format_txt = f"_{rois_method_format}"
    
            # define directories and fn
            ncsf_dir = "{}/{}/derivatives/pp_data/{}/{}/ncsf".format(
                main_dir, project_dir, subject, format_)
            
            if not os.path.isdir(ncsf_dir):
                print(f"[SKIP] ncsf_dir not found for format={format_}: {ncsf_dir}")
                continue
            
            fit_dir = "{}/fit".format(ncsf_dir)
            prf_deriv_dir = "{}/ncsf_derivatives".format(ncsf_dir)
            flatmaps_dir = '{}/pycortex/flatmaps_ncsf'.format(ncsf_dir)
            datasets_dir = '{}/pycortex/datasets_ncsf'.format(ncsf_dir)
            os.makedirs(flatmaps_dir, exist_ok=True)
            os.makedirs(datasets_dir, exist_ok=True)
        
            # get run numbers
            if "single-run" in avg_method:
                if format_ == 'fsnative':
                    deriv_list = glob.glob('{}/{}_task-{}*_hemi-L_{}_{}_{}_{}_prf-gauss_deriv.func.gii'.format(
                            prf_deriv_dir, subject, ncsf_task_name, 
                            preproc_prep, filtering, normalization, avg_method))
                elif format_ == '170k':
                    deriv_list = glob.glob('{}/{}_task-{}*_{}_{}_{}_{}_prf-gauss_deriv.dtseries.nii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, 
                        preproc_prep, filtering, normalization, avg_method))
                runs = sorted(set('_run-' + f.split('run-')[1].split('_')[0] for f in deriv_list))
            else:
                runs = ['']
                
            for run in runs:
                if format_ == 'fsnative':
                    pycortex_subject = subject
                    
                    deriv_avg_fn_L = '{}/{}_task-{}{}_hemi-L_{}_{}_{}_{}_ncsf_deriv.func.gii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method)
                    deriv_avg_fn_R = '{}/{}_task-{}{}_hemi-R_{}_{}_{}_{}_ncsf_deriv.func.gii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method)
                    deriv_results = load_surface_pycortex(L_fn=deriv_avg_fn_L, 
                                                    R_fn=deriv_avg_fn_R)
                    deriv_mat = deriv_results['data_concat']
                    
                    stats_avg_fn_L = '{}/{}_task-{}{}_hemi-L_{}_{}_{}_{}_ncsf_stats.func.gii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method)
                    stats_avg_fn_R = '{}/{}_task-{}{}_hemi-R_{}_{}_{}_{}_ncsf_stats.func.gii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method)
                    stats_results = load_surface_pycortex(L_fn=deriv_avg_fn_L, 
                                                    R_fn=deriv_avg_fn_R)
                    stats_mat = stats_results['data_concat']
                    
                elif format_ == '170k':
                    pycortex_subject = pycortex_subject_template
                    deriv_avg_fn = '{}/{}_task-{}{}_{}_{}_{}_{}_ncsf_deriv.dtseries.nii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method)
                    results = load_surface_pycortex(brain_fn=deriv_avg_fn)
                    deriv_mat = results['data_concat']
                    
                    stats_avg_fn = '{}/{}_task-{}{}_{}_{}_{}_{}_ncsf_stats.dtseries.nii'.format(
                        prf_deriv_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method)
                    stats_results = load_surface_pycortex(brain_fn=stats_avg_fn)
                    stats_mat = stats_results['data_concat']
                
                # Combine mat
                all_deriv_mat = np.concatenate((deriv_mat, stats_mat))
                
                # threshold data
                deriv_mat_th = all_deriv_mat
                if analysis_info['stats_th'] == 0.05: stats_th_down = deriv_mat_th[corr_pvalue_5pt_idx,...] <= 0.05
                elif analysis_info['stats_th'] == 0.01: stats_th_down = deriv_mat_th[corr_pvalue_1pt_idx,...] <= 0.01
                all_th = np.array((stats_th_down)) 
                all_deriv_mat[rsq_idx2use, all_th==False]=0

                # Create flatmaps
                print('Creating flatmaps...')
                maps_names = []
    
                # r-square
                rsq_data = all_deriv_mat[rsq_idx2use,...]
                alpha = (rsq_data - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
                alpha[alpha>1]=1
                param_rsq = {'data': rsq_data, 
                             'cmap': cmap_r2, 
                             'alpha': alpha,
                             'vmin': rsq_scale[0], 
                             'vmax': rsq_scale[1], 
                             'cbar': 'discrete',
                             'cortex_type': 'VertexRGB',
                             'description': rsq_description,
                             'curv_brightness': 1,
                             'curv_contrast': 0.1, 
                             'add_roi': save_svg,
                             'cbar_label': rsq_cbar_label, 
                             'overlay_fn': overlay_fn,
                             'with_labels': True}
                maps_names.append('rsq')
                
                # SFp
                SFp_data = all_deriv_mat[SFp_idx,...]
                param_SFp = {'data': SFp_data, 
                             'cmap': cmap_SFp_CSp_SFmax, 
                             'alpha': alpha,
                             'vmin': SFp_scale[0], 
                             'vmax': SFp_scale[1], 
                             'cbar': 'discrete', 
                             'cortex_type': 'VertexRGB',
                             'description': 'nCSF SFp', 
                             'curv_brightness': 1,
                             'curv_contrast': 0.1, 
                             'add_roi': save_svg, 
                             'cbar_label': 'nCSF SFp',
                             'overlay_fn': overlay_fn,
                             'with_labels': True}
                maps_names.append('SFp')
                
                # CSp
                CSp_data = all_deriv_mat[CSp_idx,...]
                param_CSp = {'data': CSp_data, 
                              'cmap': cmap_SFp_CSp_SFmax,
                              'alpha': alpha, 
                              'vmin': CSp_scale[0],
                              'vmax': CSp_scale[1],
                              'cbar': 'discrete', 
                              'cortex_type': 'VertexRGB',
                              'description': 'nCSF CSp', 
                              'curv_brightness': 1,
                              'curv_contrast': 0.1,
                              'add_roi': False,
                              'cbar_label': 'nCSF CSp',
                              'overlay_fn': overlay_fn,
                              'with_labels': True}
                maps_names.append('CSp')
                
                # auc
                auc_data = all_deriv_mat[auc_idx,...]
                param_auc = {'data': auc_data, 
                              'cmap': cmap_auc,
                              'alpha': alpha, 
                              'vmin': auc_scale[0],
                              'vmax': auc_scale[1],
                              'cbar': 'discrete', 
                              'cortex_type': 'VertexRGB',
                              'description': 'nCSF auc', 
                              'curv_brightness': 1,
                              'curv_contrast': 0.1,
                              'add_roi': False,
                              'cbar_label': 'nCSF auc',
                              'overlay_fn': overlay_fn,
                              'with_labels': True}
                maps_names.append('auc')
                
                # normalize auc
                normalize_auc_data = all_deriv_mat[normalize_auc_idx,...]
                param_normalize_auc = {'data': normalize_auc_data, 
                              'cmap': cmap_auc,
                              'alpha': alpha, 
                              'vmin': normalize_auc_scale[0],
                              'vmax': normalize_auc_scale[1],
                              'cbar': 'discrete', 
                              'cortex_type': 'VertexRGB',
                              'description': 'nCSF normalize auc', 
                              'curv_brightness': 1,
                              'curv_contrast': 0.1,
                              'add_roi': False,
                              'cbar_label': 'nCSF normalize auc',
                              'overlay_fn': overlay_fn,
                              'with_labels': True}
                maps_names.append('normalize_auc')
                
                # SFmax
                SFmax_data = all_deriv_mat[SFmax_idx,...]
                param_SFmax = {'data': SFmax_data, 
                              'cmap': cmap_SFp_CSp_SFmax,
                              'alpha': alpha, 
                              'vmin': SFmax_scale[0],
                              'vmax': SFmax_scale[1],
                              'cbar': 'discrete', 
                              'cortex_type': 'VertexRGB',
                              'description': 'nCSF SFmax', 
                              'curv_brightness': 1,
                              'curv_contrast': 0.1,
                              'add_roi': False,
                              'cbar_label': 'nCSF SFmax',
                              'overlay_fn': overlay_fn,
                              'with_labels': True}
                maps_names.append('SFmax')
                
                # CRFsplope
                crf_exp_data = all_deriv_mat[crf_exp_idx,...]
                param_crf_exp = {'data': crf_exp_data, 
                              'cmap': cmap_SFp_CSp_SFmax,
                              'alpha': alpha, 
                              'vmin': crf_exp_scale[0],
                              'vmax': crf_exp_scale[1],
                              'cbar': 'discrete', 
                              'cortex_type': 'VertexRGB',
                              'description': 'CRF slope', 
                              'curv_brightness': 1,
                              'curv_contrast': 0.1,
                              'add_roi': False,
                              'cbar_label': 'CRF slope',
                              'overlay_fn': overlay_fn,
                              'with_labels': True}
                maps_names.append('crf_exp')
                
                # draw flatmaps
                volumes = {}
                for maps_name in maps_names:
                
                    # create flatmap
                    roi_name = 'nCSF_{}'.format(maps_name)
                    roi_param = {'subject': pycortex_subject, 
                                 'roi_name': roi_name}
                    print(roi_name)
                    exec('param_{}.update(roi_param)'.format(maps_name))
                    exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name=maps_name))
                    exec("plt.savefig('{}/{}_task-{}{}_{}_{}_{}_{}{}_nCSF-{}.pdf')".format(
                        flatmaps_dir, subject, ncsf_task_name, run,
                        preproc_prep, filtering, normalization, avg_method,
                        rois_method_format_txt, maps_name))
                    plt.close()
                
                    # save flatmap as dataset
                    exec('vol_description = param_{}["description"]'.format(maps_name))
                    exec('volume = volume_{}'.format(maps_name))
                    volumes.update({vol_description:volume})
                
                # save dataset
                dataset_file = "{}/{}_task-{}{}_{}_{}_{}_{}{}_nCSF.hdf".format(
                    datasets_dir, subject, ncsf_task_name, run,
                    preproc_prep, filtering, normalization, 
                    avg_method, rois_method_format_txt)
                if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
                dataset = cortex.Dataset(data=volumes)
                dataset.save(dataset_file)
                
                
                