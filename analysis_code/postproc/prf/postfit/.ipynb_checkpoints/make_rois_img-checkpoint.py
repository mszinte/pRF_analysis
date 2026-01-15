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
python make_rois_img.py /scratch/mszinte/data RetinoMaps sub-170k 327
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
import json
import numpy as np
import pandas as pd
import nibabel as nb
import glob

# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from cifti_utils import from_170k_to_59k, from_59k_to_170k
from surface_utils import make_surface_image, load_surface
from pycortex_utils import get_rois, set_pycortex_config_file

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
rois = analysis_info["rois"]
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
if subject == 'sub-170k': extensions = ['dtseries.nii']
else: extensions = analysis_info['extensions']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Create roi image files
for format_, extension in zip(formats, extensions): 
    print(format_)
    rois_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(
        main_dir, project_dir, subject, format_)
    os.makedirs(rois_dir, exist_ok=True)
    
    if format_ == 'fsnative':        
        # Load rois 
        roi_verts_dict_L, roi_verts_dict_R = get_rois(subject, 
                                                      return_concat_hemis=False, 
                                                      rois=rois, 
                                                      mask=True, 
                                                      atlas_name=None, 
                                                      surf_size=None)
        
        for hemi in ['hemi-L','hemi-R']:
            if hemi == 'hemi-L':roi_verts_dict = roi_verts_dict_L
            elif hemi == 'hemi-R':roi_verts_dict = roi_verts_dict_R
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
            rois_fn = '{}_{}_rois.{}'.format(subject, hemi, extension)

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
        roi_verts_dict = get_rois(subject, 
                                  return_concat_hemis=True, 
                                  rois=rois,
                                  mask=True, 
                                  atlas_name='mmp_group', 
                                  surf_size='170k')

        array_rois = np.zeros(len(next(iter(roi_verts_dict.values()))), dtype=int)  
        for i, (key, mask) in enumerate(roi_verts_dict.items(), 1):
            array_rois[mask] = i
            
        # Define filename
        rois_fn = '{}_rois.{}'.format(subject, extension)

        # Saving file
        array_rois = array_rois.reshape(1, -1)
        rois_img = make_surface_image(data=array_rois, source_img=img, maps_names=['rois'])
        nb.save(rois_img, '{}/{}'.format(rois_dir, rois_fn))
        print('Saving {}/{}'.format(rois_dir, rois_fn))
        
        # MMP rois
        #  Load csv 
        sub_170k_cortex_dir = '{}/db/sub-170k'.format(cortex_dir)

        # load npz atlas and make dataframe
        mmp_atlas_fn ='{}/surface-info/mmp_atlas.npz'.format(sub_170k_cortex_dir)
        mmp_npz = np.load(mmp_atlas_fn)
        
        #  make left hemi 
        mmp_df_lh = pd.DataFrame(mmp_npz['left'], columns = ['roi_id'])
        mmp_df_lh = mmp_df_lh.assign(hemi='L')
         
        #  make right hemi 
        mmp_df_rh = pd.DataFrame(mmp_npz['right'], columns = ['roi_id'])
        mmp_df_rh = mmp_df_rh.assign(hemi='R')
        
        #  make brain df
        mmp_df_brain = pd.concat([mmp_df_lh,mmp_df_rh], ignore_index=True )
        mmp_df_brain = mmp_df_brain.assign(roi_name=np.nan)
        
        mmp_df_brain['roi_id'] = mmp_df_brain['roi_id'].replace({0: 180})
        mmp_df_brain['roi_id_hemi'] = np.where(mmp_df_brain['hemi'] == 'R', 
                                               mmp_df_brain['roi_id'] + 180, 
                                               mmp_df_brain['roi_id'])
        
        #  Load mmp information csv 
        mmp_csv_fn = '{}/HCP-MMP1_UniqueRegionList.csv'.format(sub_170k_cortex_dir)
        mmp_csv = pd.read_csv(mmp_csv_fn)
        mmp_csv['index_col'] = (mmp_csv.index + 1).astype('int32')
        
        # make the final dataframe with the correpondamce between the code and the areas
        mmp_final_df = pd.merge(mmp_df_brain, 
                                mmp_csv, 
                                left_on='roi_id_hemi', 
                                right_on='index_col', 
                                how='left')
        mmp_final_df = mmp_final_df[['roi_id', 'roi_id_hemi', 'region','hemi']]  
        mmp_final_df.rename(columns={'region': 'roi_name'}, inplace=True)
        
        mmp_array = np.array(mmp_final_df['roi_id']).reshape(1,-1)
        
        # Get 59k mask
        results = from_170k_to_59k(img=img, 
                                   data=data, 
                                   return_concat_hemis=True, 
                                   return_59k_mask=True)

        mask_59k = results['mask_59k']
        
        # convert mmp rois in 170k 
        mmp_array_170k = from_59k_to_170k(data_59k=mmp_array,
                                          brain_mask_59k=mask_59k)
        
        # Define filename
        rois_mmp_fn = '{}_rois_mmp.{}'.format(subject, extension)

        # Saving file
        array_rois = array_rois.reshape(1, -1)
        rois_img = make_surface_image(data=mmp_array_170k, source_img=img, maps_names=['rois'])
        nb.save(rois_img, '{}/{}'.format(rois_dir, rois_mmp_fn))
        print('Saving {}/{}'.format(rois_dir, rois_fn))

# Change permission
print('Changing permission in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))