"""
-----------------------------------------------------------------------------------------
compute_gauss_derivatives.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute pRF derivatives from the pRF grid gauss fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis task name (ex. prf)
sys.argv[5]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and pRF derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_gauss_derivatives.py [main directory] [project name] 
                                       [subject num] [analysis name][group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python compute_gauss_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 prf meso_0327
python compute_gauss_derivatives.py /scratch/mszinte/data RetinoMaps template_avg prf meso_0327
python compute_gauss_derivatives.py /scratch/mszinte/data amblyo7T_prf sub-01 prf meso_0327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
edited by Sina Kling (sina.kling@outlook.de)
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
import nibabel as nb

# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from prf_utils import fit2deriv
from maths_utils import  median_subject_template
from surface_utils import make_surface_image , load_surface
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
task_names = analysis_info['analysis_task_names']
maps_names_gauss = analysis_info['maps_names_gauss']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
# [NORMALIZATION] list of normalizations, looped over below
# AVANT : normalization = analysis_info['normalization']
normalizations = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
averaging_templates = analysis_info['averaging_templates']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]

# template_avg exception
if subject != 'template_avg':
    # [NORMALIZATION] loop over normalizations (inside the if/elif: the template branch reassigns `subject`)
    for normalization in normalizations:
        for format_, extension in zip(formats, extensions):
            print(format_)
            # Define directories
            pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
            prf_fit_dir = "{}/{}/{}/{}/fit".format(pp_dir, subject, format_, analysis_name)
            prf_deriv_dir = "{}/{}/{}/{}/prf_derivatives".format(pp_dir, subject, format_, analysis_name)
            os.makedirs(prf_deriv_dir, exist_ok=True)
        
            # Get prf fit filenames
            # [NORMALIZATION] glob restricted to the current normalization (otherwise every normalization is recomputed at each pass)
            # AVANT : fit_fns = glob.glob("{}/{}/{}/{}/fit/*{}-gauss{}_fit*".format(pp_dir, subject, format_, analysis_name, analysis_name, dm_name))
            fit_fns = glob.glob("{}/{}/{}/{}/fit/*_{}_{}_{}_*{}-gauss{}_fit*".format(pp_dir, subject, format_, analysis_name, preproc_prep, filtering, normalization, analysis_name, dm_name))

            # Compute derivatives 
            for fit_fn in fit_fns:
            
                deriv_fn = fit_fn.split('/')[-1]
                deriv_fn = deriv_fn.replace(f'{analysis_name}-gauss{dm_name}_fit', f'{analysis_name}-gauss{dm_name}_deriv')
        
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
                    print('Saving {}/{}'.format(prf_deriv_dir, deriv_fn))
                    nb.save(deriv_img,'{}/{}'.format(prf_deriv_dir, deriv_fn))

    # template_avg median          
elif subject == 'template_avg':
    # [NORMALIZATION] loop over normalizations (inside the if/elif: the template branch reassigns `subject`)
    for normalization in normalizations:
        for averaging_template_name, averaging_template_format in averaging_templates.items(): 
            print('{}, Median corr across subject...'. format(averaging_template_name))
        
            for task_name in task_names:
            
                for avg_method in avg_methods:
                    print('{}'.format(avg_method))
                    if "loo" in avg_method:
                        continue  # Skip if it contains "loo"
                
                    # find all the subject prf derivatives
                    prf_deriv_fns = []
                    for subject in subjects: 
                        prf_deriv_dir = '{}/{}/derivatives/pp_data/{}/{}/{}/prf_derivatives'.format(main_dir, project_dir, subject, averaging_template_format, analysis_name)
                        # [CORRECTION #17] nom de fichier explicite (un par sujet) : le glob *{avg_method}* pouvait aussi attraper des fichiers loo-avg-N
                        # AVANT : prf_deriv_fns_subject = "{}/{}_task-{}*{}*{}-gauss{}_deriv.dtseries.nii".format(prf_deriv_dir, subject, task_name, avg_method, analysis_name, dm_name)
                        # AVANT : prf_deriv_fns.extend(glob.glob(prf_deriv_fns_subject))
                        prf_deriv_fns += ["{}/{}_task-{}_{}_{}_{}_{}_{}-gauss{}_deriv.dtseries.nii".format(
                            prf_deriv_dir, subject, task_name, preproc_prep, filtering, normalization, avg_method, analysis_name, dm_name)]
                
                    # Median across subject
                    img, data_deriv_median = median_subject_template(fns=prf_deriv_fns)
    
                    # Export results
                    template_deriv_dir = "{}/{}/derivatives/pp_data/{}/{}/{}/prf_derivatives".format(
                            main_dir, project_dir, averaging_template_name, averaging_template_format, analysis_name)
                    os.makedirs(template_deriv_dir, exist_ok=True)
                
                    template_deriv_fn = "{}/{}_task-{}_{}_{}_{}_{}_{}-gauss{}_deriv.dtseries.nii".format(template_deriv_dir, averaging_template_name, task_name, preproc_prep, filtering, normalization, avg_method, analysis_name, dm_name)
                    print("saving: {}".format(template_deriv_fn))
                    template_deriv_img = make_surface_image(
                        data=data_deriv_median, source_img=img, maps_names=maps_names_gauss)
                    nb.save(template_deriv_img, template_deriv_fn)
        
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))