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
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_derivatives.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps sub-170k 327
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
from maths_utils import  avg_subject_template
from surface_utils import make_surface_image , load_surface
from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois = analysis_info['rois']
maps_names = analysis_info['maps_names_css']
subjects = analysis_info['subjects']
prf_task_name = analysis_info['prf_task_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define folders
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# sub-170k exeption
if subject != 'sub-170k':
    for format_, extension in zip(formats, extensions):
        print(format_)
        
        # Define directories
        prf_deriv_dir = "{}/{}/{}/prf/prf_derivatives".format(pp_dir, subject, format_)
        os.makedirs(prf_deriv_dir, exist_ok=True)
        
        # Get prf fit filenames
        fit_fns = glob.glob("{}/{}/{}/prf/fit/*prf-fit_css*.{}".format(
            pp_dir, subject, format_, extension))
        
        # Compute derivatives
        for fit_fn in fit_fns:
            deriv_fn = fit_fn.split('/')[-1].replace('prf-fit', 'prf-deriv')
        
            if os.path.isfile(fit_fn) == False:
                sys.exit('Missing files, analysis stopped : {}'.format(fit_fn))
            else:
                # get arrays
                fit_img, fit_data = load_surface(fit_fn)
                deriv_array = fit2deriv(fit_array=fit_data, model='css', is_loo_r2=True)
                deriv_img = make_surface_image(data=deriv_array, 
                                               source_img=fit_img, 
                                               maps_names=maps_names)
    
            nb.save(deriv_img,'{}/{}'.format(prf_deriv_dir, deriv_fn))
            print('Saving derivatives: {}'.format('{}/{}'.format(prf_deriv_dir, deriv_fn)))
    
    # Find all the derivatives files 
    derives_fns = []
    for format_, extension in zip(formats, extensions):
        list_ = glob.glob("{}/{}/{}/prf/prf_derivatives/*loo-*_prf-deriv_css.{}".format(
            pp_dir, subject, format_, extension))
        derives_fns.extend(list_)
    
    # Split filtered files depending of their nature
    deriv_fsnative_hemi_L, deriv_fsnative_hemi_R, deriv_170k = [], [], []
    for subtype in derives_fns:
        if "hemi-L" in subtype: deriv_fsnative_hemi_L.append(subtype)
        elif "hemi-R" in subtype: deriv_fsnative_hemi_R.append(subtype)
        else : deriv_170k.append(subtype)
    
    loo_deriv_fns_list = [deriv_fsnative_hemi_L,
                          deriv_fsnative_hemi_R, 
                          deriv_170k]
    hemi_data_avg = {'hemi-L': [], 
                     'hemi-R': [], 
                     '170k': []}
    
    # Averaging
    for loo_deriv_fns in loo_deriv_fns_list:
        if loo_deriv_fns[0].find('hemi-L') != -1: hemi = 'hemi-L'
        elif loo_deriv_fns[0].find('hemi-R') != -1: hemi = 'hemi-R'
        else: hemi = None
    
        deriv_img, deriv_data = load_surface(fn=loo_deriv_fns[0])
        loo_deriv_data_avg = np.zeros(deriv_data.shape)
        for n_run, loo_deriv_fn in enumerate(loo_deriv_fns):
            loo_deriv_avg_fn = loo_deriv_fn.split('/')[-1]
            loo_deriv_avg_fn = re.sub(r'avg_loo-\d+_prf-deriv', 'prf-deriv-loo-avg', loo_deriv_avg_fn)
            
            # Load data 
            loo_deriv_img, loo_deriv_data = load_surface(fn=loo_deriv_fn)
            
            # Averaging
            if n_run == 0: loo_deriv_data_avg = np.copy(loo_deriv_data)
            else: loo_deriv_data_avg = np.nanmean(np.array([loo_deriv_data_avg, loo_deriv_data]), axis=0)
        
        if hemi:
            avg_fn = '{}/{}/fsnative/prf/prf_derivatives/{}'.format(
                pp_dir, subject, loo_deriv_avg_fn)
            hemi_data_avg[hemi] = loo_deriv_data_avg
        else:
            avg_fn = '{}/{}/170k/prf/prf_derivatives/{}'.format(
                pp_dir, subject, loo_deriv_avg_fn)
            hemi_data_avg['170k'] = loo_deriv_data_avg
            
        # Export averaged data in surface format 
        loo_deriv_img = make_surface_image(data=loo_deriv_data_avg, 
                                           source_img=loo_deriv_img, 
                                           maps_names=maps_names)
        nb.save(loo_deriv_img, avg_fn)
        print('Saving avg: {}'.format(avg_fn))

# Sub-170k averaging                
elif subject == 'sub-170k':
    print('sub-170, averaging prf deriv across subject...')
    # find all the subject prf derivatives
    subjects_derivatives = []
    for subject in subjects: 
        subjects_derivatives += ["{}/{}/derivatives/pp_data/{}/170k/prf/prf_derivatives/{}_task-{}_fmriprep_dct_prf-deriv-loo-avg_css.dtseries.nii".format(
                main_dir, project_dir, subject, subject, prf_task_name)]

    # Averaging across subject
    img, data_deriv_avg = avg_subject_template(fns=subjects_derivatives)
        
    # Export results
    sub_170k_deriv_dir = "{}/{}/derivatives/pp_data/sub-170k/170k/prf/prf_derivatives".format(
            main_dir, project_dir)
    os.makedirs(sub_170k_deriv_dir, exist_ok=True)
    
    sub_170k_deriv_fn = "{}/sub-170k_task-{}_fmriprep_dct_prf-deriv-loo-avg_css.dtseries.nii".format(sub_170k_deriv_dir, prf_task_name)
    
    print("save: {}".format(sub_170k_deriv_fn))
    sub_170k_deriv_img = make_surface_image(
        data=data_deriv_avg, source_img=img, maps_names=maps_names)
    nb.save(sub_170k_deriv_img, sub_170k_deriv_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))