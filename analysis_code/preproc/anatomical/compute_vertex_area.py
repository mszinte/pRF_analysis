"""
-----------------------------------------------------------------------------------------
compute_vertex_area.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute vertex area
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# Image of vertex area
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/preproc/anatomical/
2. run python command
python compute_vertex_area.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical/

python compute_vertex_area.py /scratch/mszinte/data MotConf sub-01 327
python compute_vertex_area.py /scratch/mszinte/data MotConf sub-170k 327

python compute_vertex_area.py /scratch/mszinte/data RetinoMaps sub-01 327
python compute_vertex_area.py /scratch/mszinte/data RetinoMaps sub-170k 327

python compute_vertex_area.py /scratch/mszinte/data amblyo_prf sub-01 327
python compute_vertex_area.py /scratch/mszinte/data amblyo_prf sub-170k 327

python compute_vertex_area.py /scratch/mszinte/data centbids sub-2100247523 327
python compute_vertex_area.py /scratch/mszinte/data centbids sub-170k 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
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
import glob
import json
import cortex
import nibabel as nb

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from cifti_utils import from_59k_to_170k
from surface_utils import load_surface , make_surface_image
from pycortex_utils import calculate_vertex_areas, load_surface_pycortex, set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
maps_names_vert_area = analysis_info["maps_names_vert_area"]
subjects = analysis_info["subjects"]

# Find all the filtered files 
preproc_fns = []
for format_, extension in zip(formats, extensions):
    if subject == 'sub-170k': # pick first subject sub-170k
        list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct/*_*.{}".format(
                main_dir, project_dir, subjects[0], format_, extension))
    else:    
        list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct/*_*.{}".format(
                main_dir, project_dir, subject, format_, extension))
    preproc_fns.extend(list_)

# Split filtered files  depending of their nature
preproc_fsnative_hemi_L, preproc_fsnative_hemi_R, preproc_170k = [], [], []
for subtype in preproc_fns:
    if "hemi-L" in subtype:
        preproc_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype:
        preproc_fsnative_hemi_R.append(subtype)
    elif "170k" in subtype:
        preproc_170k.append(subtype)
        
preproc_files_list = [preproc_fsnative_hemi_L,
                      preproc_fsnative_hemi_R,
                      preproc_170k]

# compute vertex area 
for format_, extension in zip(formats, extensions): 
    print("Computing vertex area {}".format(format_))
    
    # Define output folders
    dest_dir = "{}/{}/derivatives/pp_data/{}/{}/vertex_area".format(
        main_dir, project_dir, subject, format_)
    os.makedirs(dest_dir, exist_ok=True)

    if format_ == 'fsnative':
        # Get pycortex surface
        surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(subject, "flat")]
        surf_lh, surf_rh = surfs[0], surfs[1]
        
        for hemi, surf in zip(['hemi-L', 'hemi-R'],[surf_lh,surf_rh]):
            if hemi == 'hemi-L': preproc_fn = preproc_fsnative_hemi_L[0]
            elif hemi == 'hemi-R': preproc_fn = preproc_fsnative_hemi_R[0]
            
            # Load data to get source img 
            img, data = load_surface(fn=preproc_fn)
            
            # Compute vertex area 
            vertex_area = calculate_vertex_areas(pts=surf.pts, polys=surf.polys)
            vertex_area = vertex_area.reshape(1,-1)
            
            # Save image
            vertex_area_fn = '{}/{}_{}_vertex_area.{}'.format(dest_dir, subject, hemi, extension)
            vert_surf_img = make_surface_image(data=vertex_area, 
                                               source_img=img, 
                                               maps_names=maps_names_vert_area)
            print("Saving: {}".format(vertex_area_fn))
            nb.save(vert_surf_img, vertex_area_fn)
                
    elif format_ == '170k': 
        # Load data to get source img and 59k mask 
        preproc_fn = preproc_170k[0]
        
        # Acces to a 59k mask 
        results = load_surface_pycortex(brain_fn=preproc_fn, 
                                        return_img=True, 
                                        return_59k_mask=True)
        img, mask_59k = results['img'], results['mask_59k']

        # Get pycortex surface polys and pts
        pts, polys = cortex.db.get_surf('sub-170k', "flat",  merge=True)
        
        # Compute vertex area 
        vertex_area_59k = calculate_vertex_areas(pts=pts, polys=polys)
        vertex_area_59k = vertex_area_59k.reshape(1,-1)
        
        # Converte vertex area from 59k to 170k
        vertex_area_170k = from_59k_to_170k(data_59k=vertex_area_59k, 
                                            brain_mask_59k=mask_59k)
        
        # Save image
        vertex_area_fn = '{}/{}_vertex_area.{}'.format(dest_dir, subject, extension)
        vertex_area_170k_img = make_surface_image(data=vertex_area_170k, 
                                                  source_img=img, 
                                                  maps_names=maps_names_vert_area)
        print("Saving: {}".format(vertex_area_fn))
        nb.save(vertex_area_170k_img, vertex_area_fn)
        
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))