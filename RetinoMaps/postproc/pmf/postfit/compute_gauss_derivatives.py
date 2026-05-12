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
sys.argv[4]: group (e.g. 327)
sys.argv[5]: fit type ('residuals', 'bold', or 'prf')
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and pRF derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
>> python compute_gauss_derivatives.py [main directory] [project name] 
                                       [subject num] [group] [fit_type]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python compute_gauss_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327 residuals
python compute_gauss_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327 bold
python compute_gauss_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327 prf
python compute_gauss_derivatives.py /scratch/mszinte/data RetinoMaps template_avg 327 residuals
python compute_gauss_derivatives.py /scratch/mszinte/data amblyo7T_prf sub-01 327 bold
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

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from prf_utils import fit2deriv
from maths_utils import median_subject_template
from surface_utils import make_surface_image, load_surface
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
fit_type = sys.argv[5]  # 'residuals', 'bold', or 'prf'

# Validate fit_type
valid_fit_types = ['residuals', 'bold', 'prf']
if fit_type not in valid_fit_types:
    sys.exit("Invalid fit_type '{}'. Must be one of: {}".format(fit_type, valid_fit_types))

# Define fit_type-specific filename patterns
fit_type_config = {
    'residuals': {
        'fit_pattern':   'pmf-residuals-gauss_fit',
        'deriv_tag':     'pmf-residuals-gauss_deriv',
    },
    'bold': {
        'fit_pattern':   'pmf2-gauss_fit',
        'deriv_tag':     'pmf2-gauss_deriv',
    },
    'prf': {
        'fit_pattern':   'prf-gauss_fit',
        'deriv_tag':     'prf-gauss_deriv',
    },
}

cfg = fit_type_config[fit_type]
fit_pattern = cfg['fit_pattern']
deriv_tag   = cfg['deriv_tag']

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
prf_task_names = analysis_info['prf_task_names']
maps_names_gauss = analysis_info['maps_names_gauss']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = ['concat']
averaging_templates = analysis_info['averaging_templates']

# -----------------------------------------------------------------------------------------
# Per-subject processing
# -----------------------------------------------------------------------------------------
if subject != 'template_avg':
    for format_, extension in zip(formats, extensions):
        print(format_)

        # Define directories
        pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
        prf_fit_dir  = "{}/{}/{}/pmf/fit".format(pp_dir, subject, format_)
        prf_deriv_dir = "{}/{}/{}/pmf/pmf_derivatives".format(pp_dir, subject, format_)
        os.makedirs(prf_deriv_dir, exist_ok=True)

        # Find fit files matching the selected fit_type pattern
        fit_fns = glob.glob("{}/{}/{}/pmf/fit/*{}*".format(
            pp_dir, subject, format_, fit_pattern))

        # Compute derivatives
        for fit_fn in fit_fns:
            deriv_fn = fit_fn.split('/')[-1]
            deriv_fn = deriv_fn.replace(fit_pattern, deriv_tag)

            if not os.path.isfile(fit_fn):
                sys.exit('Missing files, analysis stopped : {}'.format(fit_fn))
            else:
                print('Computing derivatives: {}'.format(deriv_fn))

                # Load fit data
                fit_img, fit_data = load_surface(fit_fn)

                # Compute and save derivative
                deriv_array = fit2deriv(fit_array=fit_data, model='gauss')
                deriv_img = make_surface_image(data=deriv_array,
                                               source_img=fit_img,
                                               maps_names=maps_names_gauss)
                nb.save(deriv_img, '{}/{}'.format(prf_deriv_dir, deriv_fn))

# -----------------------------------------------------------------------------------------
# template_avg: median across subjects
# -----------------------------------------------------------------------------------------
elif subject == 'template_avg':
    for averaging_template_name, averaging_template_format in averaging_templates.items():
        print('{}, Median corr across subjects...'.format(averaging_template_name))

        for prf_task_name in prf_task_names:
            for avg_method in avg_methods:
                if "loo" in avg_method:
                    continue

                prf_deriv_fns = []
                for subj in subjects:
                    prf_deriv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/pmf_derivatives'.format(
                        main_dir, project_dir, subj, averaging_template_format)
                    prf_deriv_fns_subject = "{}/{}_task-{}*{}*{}.dtseries.nii".format(
                        prf_deriv_dir, subj, prf_task_name, avg_method, deriv_tag)
                    prf_deriv_fns.extend(glob.glob(prf_deriv_fns_subject))

                # Median across subjects
                img, data_deriv_median = median_subject_template(fns=prf_deriv_fns)

                # Export results
                template_deriv_dir = "{}/{}/derivatives/pp_data/{}/{}/pmf/pmf_derivatives".format(
                    main_dir, project_dir, averaging_template_name, averaging_template_format)
                os.makedirs(template_deriv_dir, exist_ok=True)

                template_deriv_fn = "{}/{}_task-{}_{}_{}_{}_{}_{}_{}.dtseries.nii".format(
                    template_deriv_dir,
                    averaging_template_name,
                    prf_task_name,
                    preproc_prep,
                    filtering,
                    normalization,
                    avg_method,
                    fit_type,       # e.g. 'residuals', 'bold', 'prf'
                    deriv_tag)

                print("Saving: {}".format(template_deriv_fn))
                template_deriv_img = make_surface_image(
                    data=data_deriv_median,
                    source_img=img,
                    maps_names=maps_names_gauss)
                nb.save(template_deriv_img, template_deriv_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(
    main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(
    main_dir=main_dir, project_dir=project_dir, group=group))