"""
-----------------------------------------------------------------------------------------
compute_css_pcm.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute population cortical magnification and add to derivatives
Note: 
CM is computed using the geodesic distances (mm) of vertices located within a radius on
the flatten surface (see vertex_cm_rad) and restricted by the ROI boundaries
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
New brain volume in derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit
2. run python command
>> python compute_pcm.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
python compute_css_pcm.py /scratch/mszinte/data MotConf sub-01 327
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
import glob
import cortex
import datetime
import numpy as np
import nibabel as nb

# Personal iports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import make_surface_image 
from maths_utils import avg_subject_template, weighted_nan_mean, weighted_nan_median
from pycortex_utils import set_pycortex_config_file, load_surface_pycortex, get_rois, make_image_pycortex

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Define analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

task = analysis_info["prf_task_name"]
vert_dist_th = analysis_info['vertex_pcm_rad']
formats = analysis_info['formats']
rois = analysis_info["rois"]
maps_names = analysis_info['maps_names_pcm']
subjects = analysis_info['subjects']
prf_task_name = analysis_info['prf_task_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Derivatives idx 
rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx, size_idx = 0, 1, 2, 3, 4
amp_idx, baseline_idx, x_idx, y_idx, hrf_1_idx = 5, 6, 7, 8, 9
hrf_2_idx, n_idx, loo_rsq_idx = 10, 11, 12

# Stats idx
slope_idx, intercept_idx, rvalue_idx, pvalue_idx = 0, 1, 2, 3
stderr_idx, trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 4, 5, 6, 7

# compute duration 
start_time = datetime.datetime.now()

# sub-170k exeption
if subject != 'sub-170k':
    # Compute PCM
    for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
        print(format_)
        
        # define directories and fn
        prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, subject, format_)
        fit_dir = "{}/fit".format(prf_dir)
        prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    
        if format_ == 'fsnative':
            # Derivatives
            atlas_name = None 
            surf_size = None        
            deriv_avg_fn_L = glob.glob('{}/{}*hemi-L*prf-deriv-*avg_css.func.gii'.format(
                prf_deriv_dir, subject))
            
            deriv_avg_fn_R = glob.glob('{}/{}*hemi-R*prf-deriv-*avg_css.func.gii'.format(
                prf_deriv_dir, subject))
               
            results = load_surface_pycortex(L_fn=deriv_avg_fn_L[0], 
                                            R_fn=deriv_avg_fn_R[0], 
                                            return_img=True)
            
            deriv_mat = results['data_concat'] 
            img_L = results['img_L'] 
            img_R = results['img_R']
            
            # Stats
            stats_avg_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name)
            stats_avg_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name)
            stats_results = load_surface_pycortex(L_fn=stats_avg_fn_L, R_fn=stats_avg_fn_R)
            stats_mat = stats_results['data_concat']
            
        elif format_ == '170k':
            # Derivatives
            atlas_name = 'mmp_group'
            surf_size = '59k'
            deriv_avg_fn = glob.glob('{}/{}*prf-deriv-*avg_css.dtseries.nii'.format(
                prf_deriv_dir, subject))
            
            results = load_surface_pycortex(brain_fn=deriv_avg_fn[0],
                                            return_img=True,
                                            return_59k_mask=True,  
                                            return_source_data=True)
            
            deriv_mat = results['data_concat']
            mask_59k = results['mask_59k']
            deriv_mat_170k = results['source_data'] 
            img = results['img']
            
            # Stats
            stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_loo-avg_prf-stats.dtseries.nii'.format(
                prf_deriv_dir, subject, prf_task_name)
            stats_results = load_surface_pycortex(brain_fn=stats_avg_fn)
            stats_mat = stats_results['data_concat']
        
        # Threshold data
        deriv_mat_th = deriv_mat
        amp_down = deriv_mat_th[amp_idx,...] > 0
        rsq_down = deriv_mat_th[loo_rsq_idx,...] >= analysis_info['rsqr_th']
        size_th_down = deriv_mat_th[size_idx,...] >= analysis_info['size_th'][0]
        size_th_up = deriv_mat_th[size_idx,...] <= analysis_info['size_th'][1]
        ecc_th_down = deriv_mat_th[ecc_idx,...] >= analysis_info['ecc_th'][0]
        ecc_th_up = deriv_mat_th[ecc_idx,...] <= analysis_info['ecc_th'][1]
        n_th_down = deriv_mat_th[n_idx,...] >= analysis_info['n_th'][0]
        n_th_up = deriv_mat_th[n_idx,...] <= analysis_info['n_th'][1]
        if analysis_info['stats_th'] == 0.05: stats_th_down = deriv_mat_th[corr_pvalue_5pt_idx,...] <= 0.05
        elif analysis_info['stats_th'] == 0.01: stats_th_down = deriv_mat_th[corr_pvalue_1pt_idx,...] <= 0.01
        all_th = np.array((amp_down,
                           rsq_down,
                           size_th_down,size_th_up, 
                           ecc_th_down, ecc_th_up,
                           n_th_down, n_th_up,
                           stats_th_down)) 
        th_idx = np.where(np.logical_and.reduce(all_th)==False)[0]

        # Get surfaces for each hemisphere
        surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(pycortex_subject, "flat")]
        surf_lh, surf_rh = surfs[0], surfs[1]
        
        # Get the vertices number per hemisphere
        lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]
        vert_num = lh_vert_num + rh_vert_num
        
        # Get a dict with the surface vertices contained in each ROI
        roi_verts_dict = get_rois(pycortex_subject, 
                                  return_concat_hemis=True, 
                                  rois=rois,
                                  mask=False, 
                                  atlas_name=atlas_name, 
                                  surf_size=surf_size)
        
        # Exclude vertex out of threshold
        th_idx_set = set(th_idx)
        for roi, indices in roi_verts_dict.items():
            # Filter out the indices that are in th_idx_set
            filtered_indices = [idx for idx in indices if idx not in th_idx_set]
            # Update the dictionary with the filtered indices
            roi_verts_dict[roi] = np.array(filtered_indices)
        
        # Derivatives settings
        vert_rsq_data = deriv_mat[loo_rsq_idx, ...]
        vert_x_data = deriv_mat[x_idx, ...]
        vert_y_data = deriv_mat[y_idx, ...]
        vert_size_data = deriv_mat[size_idx, ...]
        vert_ecc_data = deriv_mat[ecc_idx, ...]
        
        # Create empty results
        vert_cm = np.zeros((7,vert_num))*np.nan
        for roi in rois:
            # Find ROI vertex
            roi_vert_lh_idx = roi_verts_dict[roi][roi_verts_dict[roi] < lh_vert_num]
            roi_vert_rh_idx = roi_verts_dict[roi][roi_verts_dict[roi] >= lh_vert_num]
            roi_surf_lh_idx = roi_vert_lh_idx
            roi_surf_rh_idx = roi_vert_rh_idx - lh_vert_num
        
            # Get median distance of surounding vertices included in threshold
            vert_lh_rsq, vert_lh_size = vert_rsq_data[:lh_vert_num], vert_size_data[:lh_vert_num]
            vert_lh_x, vert_lh_y = vert_x_data[:lh_vert_num], vert_y_data[:lh_vert_num]
            vert_rh_rsq, vert_rh_size = vert_rsq_data[lh_vert_num:], vert_size_data[lh_vert_num:]
            vert_rh_x, vert_rh_y = vert_x_data[lh_vert_num:], vert_y_data[lh_vert_num:]
        
            for hemi in ['lh', 'rh']:
                if hemi == 'lh':
                    surf = surf_lh
                    roi_vert_idx, roi_surf_idx = roi_vert_lh_idx, roi_surf_lh_idx
                    vert_rsq, vert_x, vert_y, vert_size = vert_lh_rsq, vert_lh_x, vert_lh_y, vert_lh_size
                elif hemi == 'rh':
                    surf = surf_rh
                    roi_vert_idx, roi_surf_idx = roi_vert_rh_idx, roi_surf_rh_idx
                    vert_rsq, vert_x, vert_y, vert_size = vert_rh_rsq, vert_rh_x, vert_rh_y, vert_rh_size
        
                print('ROI -> {} / Hemisphere -> {}'.format(roi, hemi))
                for i, (vert_idx, surf_idx) in enumerate(zip(roi_vert_idx, roi_surf_idx)):
                    if vert_rsq[surf_idx] > 0:
                        # Get geodesic distances (mm)
                        try :
                            geo_patch = surf.get_geodesic_patch(radius=vert_dist_th, 
                                                                vertex=surf_idx)
                        except Exception as e:
                            print("Vertex #{}: error: {} within {} mm".format(vert_idx, e, vert_dist_th))
                            geo_patch['vertex_mask'] = np.zeros(surf.pts.shape[0]).astype(bool)
                            geo_patch['geodesic_distance'] = []
        
                        vert_dist_th_idx  = geo_patch['vertex_mask']
                        vert_dist_th_dist = np.ones_like(vert_dist_th_idx)*np.nan
                        vert_dist_th_dist[vert_dist_th_idx] = geo_patch['geodesic_distance']
        
                        # Exclude vextex out of roi
                        vert_dist_th_not_in_roi_idx = [idx for idx in np.where(vert_dist_th_idx)[0] if idx not in roi_surf_idx]
                        vert_dist_th_idx[vert_dist_th_not_in_roi_idx] = False
                        vert_dist_th_dist[vert_dist_th_not_in_roi_idx] = np.nan
                        
                        # Compute the n neighbor of each vertex
                        n_neighbor = np.where(vert_dist_th_idx)[0].shape[0]
                        vert_cm[1,vert_idx] = n_neighbor

                        if np.sum(vert_dist_th_idx) > 1:

                            # Get prf parameters of vertices in geodesic distance threshold
                            vert_ctr_x, vert_ctr_y = vert_x[surf_idx], vert_y[surf_idx]
                            vert_dist_th_idx[surf_idx] = False
                            
                            # median
                            # Compute median geodesic distance 
                            vert_geo_dist_median = weighted_nan_median(vert_dist_th_dist[vert_dist_th_idx], vert_rsq[vert_dist_th_idx])
                            vert_srd_x = weighted_nan_median(vert_x[vert_dist_th_idx], vert_rsq[vert_dist_th_idx])
                            vert_srd_y = weighted_nan_median(vert_y[vert_dist_th_idx], vert_rsq[vert_dist_th_idx])
                            
                            # Compute prf center suround distance (deg)
                            vert_prf_dist = np.sqrt((vert_ctr_x - vert_srd_x)**2 + (vert_ctr_y - vert_srd_y)**2)
                            
                            # Compute cortical magnification in mm/deg (surface distance / pRF positon distance)
                            vert_cm[0,vert_idx] = vert_geo_dist_median/vert_prf_dist
                            
                            # export median geodesic and prf distance
                            vert_cm[2,vert_idx] = vert_geo_dist_median
                            vert_cm[3,vert_idx] = vert_prf_dist
                            
                            # mean
                            # Compute median geodesic distance 
                            vert_geo_dist_mean = weighted_nan_mean(vert_dist_th_dist[vert_dist_th_idx], vert_rsq[vert_dist_th_idx])

                            vert_srd_x = weighted_nan_mean(vert_x[vert_dist_th_idx], vert_rsq[vert_dist_th_idx])
                            vert_srd_y = weighted_nan_mean(vert_y[vert_dist_th_idx], vert_rsq[vert_dist_th_idx])
                            
                            # Compute prf center suround distance (deg)
                            vert_prf_dist = np.sqrt((vert_ctr_x - vert_srd_x)**2 + (vert_ctr_y - vert_srd_y)**2)
                            
                            # Compute cortical magnification in mm/deg (surface distance / pRF positon distance)
                            vert_cm[4,vert_idx] = vert_geo_dist_mean/vert_prf_dist
                            
                            # export median geodesic and prf distance
                            vert_cm[5,vert_idx] = vert_geo_dist_mean
                            vert_cm[6,vert_idx] = vert_prf_dist


        deriv_mat_new = np.zeros((4, deriv_mat.shape[1])) * np.nan
        deriv_mat_new = vert_cm
        
        # Save data
        if format_ == 'fsnative':
            # Save as image
            new_img_L, new_img_R  = make_image_pycortex(data=deriv_mat_new, 
                                                        maps_names=maps_names,
                                                        img_L=img_L, 
                                                        img_R=img_R, 
                                                        lh_vert_num=lh_vert_num, 
                                                        rh_vert_num=rh_vert_num, 
                                                        img=None, 
                                                        brain_mask_59k=None)
    
            deriv_avg_fn_L = deriv_avg_fn_L[0].split('/')[-1].replace('deriv', 'pcm')
            deriv_avg_fn_R = deriv_avg_fn_R[0].split('/')[-1].replace('deriv', 'pcm')
            print('Saving {}'.format(deriv_avg_fn_L))
            nb.save(new_img_L, '{}/{}'.format(prf_deriv_dir, deriv_avg_fn_L))
            print('Saving {}'.format(deriv_avg_fn_R))
            nb.save(new_img_R, '{}/{}'.format(prf_deriv_dir, deriv_avg_fn_R))
                    
        elif format_ == '170k':
            # Save as image
            new_img = make_image_pycortex(data=deriv_mat_new, 
                                  maps_names=None,
                                  img_L=None, 
                                  img_R=None, 
                                  lh_vert_num=None, 
                                  rh_vert_num=None, 
                                  img=img, 
                                  brain_mask_59k=mask_59k)
    
            deriv_avg_fn = deriv_avg_fn[0].split('/')[-1].replace('deriv', 'pcm')
            print('Saving {}'.format(deriv_avg_fn))
            nb.save(new_img, '{}/{}'.format(prf_deriv_dir, deriv_avg_fn))
            
# Sub-170k averaging                
elif subject == 'sub-170k':
    print('sub-170, averaging pcm across subject...')
    # find all the subject prf derivatives
    subjects_pcm = []
    for subject in subjects: 
        subjects_pcm += ["{}/{}/derivatives/pp_data/{}/170k/prf/prf_derivatives/{}_task-{}_fmriprep_dct_prf-pcm-loo-avg_css.dtseries.nii".format(
                main_dir, project_dir, subject, subject, prf_task_name)]

    # Averaging across subject
    img, data_pcm_avg = avg_subject_template(fns=subjects_pcm)
        
    # Export results
    sub_170k_pcm_dir = "{}/{}/derivatives/pp_data/sub-170k/170k/prf/prf_derivatives".format(
            main_dir, project_dir)
    os.makedirs(sub_170k_pcm_dir, exist_ok=True)
    
    sub_170k_pcm_fn = "{}/sub-170k_task-{}_fmriprep_dct_prf-pcm-loo-avg_css.dtseries.nii".format(sub_170k_pcm_dir, prf_task_name)
    
    print("save: {}".format(sub_170k_pcm_fn))
    
    sub_170k_pcm_img = make_surface_image(
        data=data_pcm_avg, source_img=img, maps_names=maps_names)
    nb.save(sub_170k_pcm_img, sub_170k_pcm_fn)

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(start_time=start_time, 
                                                                                    end_time=end_time, 
                                                                                    dur=end_time - start_time))

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))