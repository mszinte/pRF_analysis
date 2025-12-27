"""
-----------------------------------------------------------------------------------------
compute_css_derivatives.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute pRF derivatives from the pRF grid gauss fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and pRF derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_derivatives.py [main directory] [project name] 
                                     [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
import re
import sys
import glob
import json
import numpy as np
import nibabel as nb

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from prf_utils import fit2deriv
from maths_utils import  median_subject_template
from surface_utils import make_surface_image , load_surface
from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois = analysis_info['rois']
maps_names = analysis_info['maps_names_css']
subjects = analysis_info['subjects']
prf_task_names = analysis_info['prf_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define folders
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# sub-170k exception
if subject != 'sub-170k':
    for avg_method in avg_methods:
        for format_, extension in zip(formats, extensions):

            # define/create folders
            prf_fit_dir = '{}/{}/{}/prf/fit'.format(
                    pp_dir, subject, format_)
            prf_deriv_dir = "{}/{}/{}/prf/prf_derivatives".format(
                pp_dir, subject, format_)

            for prf_task_name in prf_task_names:
                print(f'{avg_method} - {format_} - {prf_task_name}')

                # Find pRF func/pred files
                fit_fns = glob.glob('{}/*task-{}*_{}*_prf-css_fit.{}'.format(
                        prf_fit_dir, prf_task_name, avg_method, extension))
    
                # Compute derivatives
                for fit_fn in fit_fns:
                    deriv_fn = fit_fn.split('/')[-1].replace('prf-css_fit', 'prf-css_deriv')

                    # get arrays
                    fit_img, fit_data = load_surface(fit_fn)
                    deriv_array = fit2deriv(fit_array=fit_data, 
                                            model='css', 
                                            is_loo_r2=True)
                    deriv_img = make_surface_image(data=deriv_array, 
                                                   source_img=fit_img, 
                                                   maps_names=maps_names)

                    nb.save(deriv_img,'{}/{}'.format(prf_deriv_dir, deriv_fn))
                    print('Saving derivatives: {}'.format('{}/{}'.format(prf_deriv_dir, deriv_fn)))
            
                # Compute median across leave-one-out fit
                if 'loo-avg' in avg_method:
                    print('Compute median across LOO')                
                    
                    # Get LOO files (excluding any with "median" in the name)
                    loo_prf_deriv_fns = glob.glob(f"{prf_deriv_dir}/*task-{prf_task_name}*_loo-avg-*_prf-css_deriv.{extension}")

                    # Group files by hemisphere/format
                    loo_prf_deriv_fsnative_hemi_L_fns = [fn for fn in loo_prf_deriv_fns if "hemi-L" in fn]
                    loo_prf_deriv_fsnative_hemi_R_fns = [fn for fn in loo_prf_deriv_fns if "hemi-R" in fn]
                    loo_prf_deriv_170k_fns = [fn for fn in loo_prf_deriv_fns if "hemi-L" not in fn and "hemi-R" not in fn]

                    # Process each group
                    for group_files, hemi in [(loo_prf_deriv_fsnative_hemi_L_fns, "_hemi-L"),
                                              (loo_prf_deriv_fsnative_hemi_R_fns, "_hemi-R"),
                                              (loo_prf_deriv_170k_fns, "")]:
                        if len(group_files)>0:
                            
                            # Load first file to initialize median array and define fn
                            prf_deriv_img, prf_deriv_data = load_surface(group_files[0])
                            loo_prf_deriv = np.zeros_like(prf_deriv_data)
                            loo_prf_deriv_fn =  '{}/{}_task-{}{}_{}_{}_{}_loo-avg_prf-css_deriv.{}'.format(
                                prf_deriv_dir, subject, prf_task_name, hemi, 
                                preproc_prep, filtering, normalization, extension)
                            
                            # Compute median across LOO runs
                            for n_run, loo_deriv_fn in enumerate(group_files):
                                print(f'Loadding loo deriv: {loo_deriv_fn}')
                                _, loo_prf_deriv_run_data = load_surface(loo_deriv_fn)
                                if n_run == 0: loo_prf_deriv = np.copy(loo_prf_deriv_run_data)
                                else: loo_prf_deriv = np.nanmedian(np.array([loo_prf_deriv, loo_prf_deriv_run_data]), axis=0)
                            
                            # Save median results
                            loo_prf_deriv_img = make_surface_image(loo_prf_deriv, prf_deriv_img, maps_names)
                            nb.save(loo_prf_deriv_img, loo_prf_deriv_fn)
                            print(f"Saving median: {loo_prf_deriv_fn}")


# # Sub-170k computing median       
# elif subject == 'sub-170k':
#     print('sub-170, computing median prf deriv across subject...')
    
#     # find all the subject prf derivatives
#     subjects_derivatives = []
#     for subject in subjects: 
#         subjects_derivatives += ["{}/{}/derivatives/pp_data/{}/170k/prf/prf_derivatives/{}_task-{}_fmriprep_dct_avg_prf-deriv_css_loo-median.dtseries.nii".format(
#                 main_dir, project_dir, subject, subject, prf_task_name)]

#     # Computing median across subject
#     img, data_deriv_median = median_subject_template(fns=subjects_derivatives)
        
#     # Export results
#     sub_170k_deriv_dir = "{}/{}/derivatives/pp_data/sub-170k/170k/prf/prf_derivatives".format(
#             main_dir, project_dir)
#     os.makedirs(sub_170k_deriv_dir, exist_ok=True)
    
#     sub_170k_deriv_fn = "{}/sub-170k_task-{}_fmriprep_dct_avg_prf-deriv_css_loo-median.dtseries.nii".format(sub_170k_deriv_dir, prf_task_name)
    
#     print("save: {}".format(sub_170k_deriv_fn))
#     sub_170k_deriv_img = make_surface_image(data=data_deriv_median, 
#                                             source_img=img, 
#                                             maps_names=maps_names)
#     nb.save(sub_170k_deriv_img, sub_170k_deriv_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))