"""
-----------------------------------------------------------------------------------------
compute_gauss_gridfit_derivatives.py
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
>> python compute_gauss_gridfit_derivatives.py [main directory] [project name] 
                                               [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python compute_gauss_gridfit_derivatives.py /scratch/mszinte/data MotConf sub-01 327
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
import json
import nibabel as nb

# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from prf_utils import fit2deriv
from maths_utils import  median_subject_template
from surface_utils import make_surface_image , load_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
prf_task_name = analysis_info['prf_task_name']
maps_names_gauss = analysis_info['maps_names_gauss']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
gauss_avg = analysis_info['gauss_avg']

# sub-170k exception
if subject != 'sub-170k':
    for format_, extension in zip(formats, extensions):
        print(format_)
        # Define directories
        pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
        prf_fit_dir = "{}/{}/{}/prf/fit".format(pp_dir, subject, format_)
        prf_deriv_dir = "{}/{}/{}/prf/prf_derivatives".format(pp_dir, subject, format_)
        os.makedirs(prf_deriv_dir, exist_ok=True)
        
        # Get prf fit filenames
        fit_fns= glob.glob("{}/{}/{}/{}/fit/*prf-gauss_fit*".format(pp_dir, subject, format_))

        # Compute derivatives 
        for fit_fn in fit_fns:
            
            deriv_fn = fit_fn.split('/')[-1]
            deriv_fn = deriv_fn.replace('prf-gauss_fit', 'prf-gauss_deriv')
        
            if os.path.isfile(fit_fn) == False:
                sys.exit('Missing files, analysis stopped : {}'.format(fit_fn))
            else:
                print('Computing derivatives: {}'.format(deriv_fn))
                
                # get arrays
                fit_img, fit_data = load_surface(fit_fn)
         
                # compute and save derivatives array
                deriv_array = fit2deriv(fit_array=fit_data, model='gauss')
                deriv_img = make_surface_image(data=deriv_array, 
                                               source_img=fit_img, 
                                               maps_names=maps_names_gauss)
                nb.save(deriv_img,'{}/{}'.format(prf_deriv_dir, deriv_fn))

# Sub-170k median                
elif subject == 'sub-170k':
    print('sub-170, computing median derivatives across subject...')
    
    # find all the subject prf derivatives
    subjects_derivatives = []
    for subject in subjects: 
        subjects_derivatives += ["{}/{}/derivatives/pp_data/{}/170k/prf/prf_derivatives/{}_task-{}_{}_{}_{}_prf-gauss_deriv.dtseries.nii".format(
                main_dir, project_dir, subject, subject, prf_task_name, preproc_prep, filtering, gauss_avg]

    # Median across subject
    img, data_deriv_median = median_subject_template(fns=subjects_derivatives)
    
    # Export results
    sub_170k_deriv_dir = "{}/{}/derivatives/pp_data/sub-170k/170k/prf/prf_derivatives".format(
            main_dir, project_dir)
    os.makedirs(sub_170k_deriv_dir, exist_ok=True)
    
    sub_170k_deriv_fn = "{}/sub-170k_task-{}_{}_{}_{}_prf-gauss_deriv.dtseries.nii".format(
        sub_170k_deriv_dir, prf_task_name, preproc_prep, filtering, gauss_avg)
    
    print("save: {}".format(sub_170k_deriv_fn))
    sub_170k_deriv_img = make_surface_image(
        data=data_deriv_median, source_img=img, maps_names=maps_names_gauss)
    nb.save(sub_170k_deriv_img, sub_170k_deriv_fn)
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))