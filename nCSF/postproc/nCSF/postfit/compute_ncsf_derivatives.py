"""
-----------------------------------------------------------------------------------------
compute_css_derivatives.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute ncsf derivatives from the ncsf grid gauss fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and ncsf derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/ncsf_analysis/analysis_code/postproc/ncsf/postfit/
2. run python command
>> python compute_css_derivatives.py [main directory] [project name] 
                                     [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/ncsf_analysis/analysis_code/postproc/ncsf/postfit/
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps template_avg 327
python compute_css_derivatives.py /scratch/mszinte/data amblyo7T_ncsf sub-01 327
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
import sys
import glob
import numpy as np
import nibabel as nb

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings
from maths_utils import  median_subject_template
from pycortex_utils import set_pycortex_config_file
from surface_utils import make_surface_image , load_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
ncsf_settings_path = os.path.join(base_dir, project_dir, "ncsf-analysis.yml")
settings = load_settings([settings_path, ncsf_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
subjects = analysis_info['subjects']
filtering = analysis_info['filtering']
extensions = analysis_info['extensions']
avg_methods = analysis_info['avg_methods']
preproc_prep = analysis_info['preproc_prep']
normalization = analysis_info['normalization']
ncsf_task_name = analysis_info['ncsf_task_name']
ncsf_maps_names = analysis_info['ncsf_maps_names']
averaging_templates = analysis_info['averaging_templates']


# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define folders
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# template_avg exception
if subject != 'template_avg':
    for avg_method in avg_methods:
        if 'loo' in avg_method:
            is_loo_r2 = True
            ncsf_maps_names = ncsf_maps_names + ['r_squared']
        else: 
            is_loo_r2 = False
            ncsf_maps_names = ncsf_maps_names
            
        for format_, extension in zip(formats, extensions):

            # define/create folders
            ncsf_fit_dir = '{}/{}/{}/ncsf/fit'.format(
                    pp_dir, subject, format_)
            ncsf_deriv_dir = "{}/{}/{}/ncsf/ncsf_derivatives".format(
                pp_dir, subject, format_)

            # Compute median across leave-one-out fit
            if 'loo-avg' in avg_method:
                print('Compute median across LOO')                
                
                # Get LOO files (excluding any with "median" in the name)
                loo_ncsf_deriv_fns = glob.glob(f"{ncsf_deriv_dir}/*task-{ncsf_task_name}_*loo-avg-*_ncsf_deriv.{extension}")

                # Group files by hemisphere/format
                loo_ncsf_deriv_fsnative_hemi_L_fns = [fn for fn in loo_ncsf_deriv_fns if "hemi-L" in fn]
                loo_ncsf_deriv_fsnative_hemi_R_fns = [fn for fn in loo_ncsf_deriv_fns if "hemi-R" in fn]
                loo_ncsf_deriv_170k_fns = [fn for fn in loo_ncsf_deriv_fns if "hemi-L" not in fn and "hemi-R" not in fn]

                # Process each group
                for group_files, hemi in [(loo_ncsf_deriv_fsnative_hemi_L_fns, "_hemi-L"),
                                          (loo_ncsf_deriv_fsnative_hemi_R_fns, "_hemi-R"),
                                          (loo_ncsf_deriv_170k_fns, "")]:
                    if len(group_files)>0:
                        
                        # Load first file to initialize median array and define fn
                        ncsf_deriv_img, ncsf_deriv_data = load_surface(group_files[0])
                        loo_ncsf_deriv = np.zeros_like(ncsf_deriv_data)
                        loo_ncsf_deriv_fn =  '{}/{}_task-{}{}_{}_{}_{}_loo-avg_ncsf-css_deriv.{}'.format(
                            ncsf_deriv_dir, subject, ncsf_task_name, hemi, 
                            preproc_prep, filtering, normalization, extension)
                        
                        # Compute median across LOO runs
                        for n_run, loo_deriv_fn in enumerate(group_files):
                            print(f'Loadding loo deriv: {loo_deriv_fn}')
                            _, loo_ncsf_deriv_run_data = load_surface(loo_deriv_fn)
                            if n_run == 0: loo_ncsf_deriv = np.copy(loo_ncsf_deriv_run_data)
                            else: loo_ncsf_deriv = np.nanmedian(np.array([loo_ncsf_deriv, loo_ncsf_deriv_run_data]), axis=0)
                        
                        # Save median results
                        loo_ncsf_deriv_img = make_surface_image(loo_ncsf_deriv, ncsf_deriv_img, ncsf_maps_names)
                        nb.save(loo_ncsf_deriv_img, loo_ncsf_deriv_fn)
                        print(f"Saving median: {loo_ncsf_deriv_fn}")


# template_avg median          
elif subject == 'template_avg':
    for averaging_template_name, averaging_template_format in averaging_templates.items(): 
        print('{}, Median corr across subject...'. format(averaging_template_name))
    
        for avg_method in avg_methods:
            if 'loo' in avg_method:
                ncsf_maps_names = ncsf_maps_names + ['r_squared']
            else: 
                ncsf_maps_names = ncsf_maps_names
            
            # find all the subject ncsf deriv
            ncsf_deriv_fns = []
            for subject in subjects:
                ncsf_deriv_dir = "{}/{}/derivatives/pp_data/{}/{}/ncsf/ncsf_derivatives".format(
                    main_dir, project_dir, subject, averaging_template_format)
                ncsf_deriv_fns += ["{}/{}_task-{}_{}_{}_{}_{}_ncsf-css_deriv.dtseries.nii".format(
                    ncsf_deriv_dir, subject, ncsf_task_name,
                    preproc_prep, filtering, normalization, avg_method)]

            # Computing median across subject
            img, data_deriv_median = median_subject_template(fns=ncsf_deriv_fns)
            
            # Export results
            template_deriv_dir = "{}/{}/derivatives/pp_data/{}/{}/ncsf/ncsf_derivatives".format(
                main_dir, project_dir, averaging_template_name, averaging_template_format)
            os.makedirs(template_deriv_dir, exist_ok=True)

            template_deriv_fn = "{}/{}_task-{}_{}_{}_{}_{}_ncsf-css_deriv.dtseries.nii".format(
                template_deriv_dir, averaging_template_name, ncsf_task_name, 
                preproc_prep, filtering, normalization, avg_method)
            
            print("save: {}".format(template_deriv_fn))
            template_deriv_img = make_surface_image(data=data_deriv_median, 
                                                    source_img=img, 
                                                    ncsf_maps_names=ncsf_maps_names)
            nb.save(template_deriv_img, template_deriv_fn)

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))