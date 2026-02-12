"""
-----------------------------------------------------------------------------------------
make_rois_img.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make Cifti and Gifti object with rois 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and pRF derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
2. run python command
>> python make_rois_img.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
python make_rois_img.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_rois_img.py /scratch/mszinte/data RetinoMaps hcp1.6mm 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import glob
import numpy as np
import nibabel as nb


# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings
from surface_utils import make_surface_image, load_surface
from pycortex_utils import get_rois, set_pycortex_config_file


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
if subject == 'sub-170k': extensions = ['dtseries.nii']
else: extensions = analysis_info['extensions']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Create roi image files
for format_, extension in zip(formats, extensions): 
    
    # define list of rois for each format
    rois_methods_format = rois_methods[format_]
    for rois_method_format in rois_methods_format:

        print(format_)
        rois_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(
            main_dir, project_dir, subject, format_)
        os.makedirs(rois_dir, exist_ok=True)

        if rois_method_format == 'rois-drawn':
            rois = analysis_info[rois_method_format]
        elif rois_method_format == 'rois-group-mmp':
            rois = list(analysis_info[rois_method_format].keys())

        if format_ == 'fsnative':                    
            for hemi in ['hemi-L','hemi-R']:
                
                roi_verts_dict = get_rois(subject=subject, 
                                          surf_format=format_, 
                                          rois_type=rois_method_format, 
                                          mask=True, 
                                          rois=rois, 
                                          hemis=hemi) 
  
                array_rois = np.zeros(len(next(iter(roi_verts_dict.values()))), dtype=int)  
                for i, (key, mask) in enumerate(roi_verts_dict.items(), 1):
                    array_rois[mask] = i
                    
                # Load data to have source img
                data_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/prf_derivatives'.format(
                    main_dir, project_dir, subject, format_)
    
                # Find first file with prf-deriv in the name
                data_files = glob.glob('{}/{}_*{}_*deriv*.{}'.format(data_dir, subject, hemi, extension))
                if not data_files:
                    raise FileNotFoundError(f"No prf-deriv file found for {subject} {hemi}")
                data_fn = data_files[0]
                img, data = load_surface(fn=data_fn)
                
                # Define filename
                rois_fn = '{}_{}_{}_{}_{}_{}.{}'.format(subject, hemi,preproc_prep, filtering, 
                                                        normalization, rois_method_format,
                                                        extension
                                                       )
    
                # Saving file
                array_rois = array_rois.reshape(1, -1)
                rois_img = make_surface_image(data=array_rois, source_img=img, maps_names=['rois'])
                nb.save(rois_img, '{}/{}'.format(rois_dir, rois_fn))
                print('Saving {}/{}'.format(rois_dir, rois_fn))
            
                
        elif format_ == '170k':
            # Load data to have source img
            data_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/prf_derivatives'.format(
                main_dir, project_dir, subject, format_)
            
            # Find first file with prf-deriv in the name
            data_files = glob.glob('{}/{}_*deriv*.{}'.format(data_dir, subject, extension))
            if not data_files:
                raise FileNotFoundError(f"No prf-deriv file found for {subject}")
            data_fn = data_files[0]
            img, data = load_surface(fn=data_fn)
            
            # MMP group rois
            roi_verts_dict = get_rois(subject=pycortex_subject_template, 
                                      surf_format=format_, 
                                      rois_type=rois_method_format,
                                      mask=True, 
                                      rois=rois, 
                                      hemis=None
                                     )

            array_rois = np.zeros(len(next(iter(roi_verts_dict.values()))), dtype=int)  
            for i, (key, mask) in enumerate(roi_verts_dict.items(), 1):
                array_rois[mask] = i
                
            # Define filename
            rois_fn = '{}_{}_{}_{}_{}.{}'.format(subject, preproc_prep, filtering, 
                                                 normalization, rois_method_format, extension)
    
            # Saving file
            array_rois = array_rois.reshape(1, -1)
            rois_img = make_surface_image(data=array_rois, source_img=img, maps_names=['rois'])
            nb.save(rois_img, '{}/{}'.format(rois_dir, rois_fn))
            print('Saving {}/{}'.format(rois_dir, rois_fn))
            
# Change permission
print('Changing permission in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))