"""
-----------------------------------------------------------------------------------------
pycortex_webgl_retinomaps.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create combined webgl per participants for pRF pycortex_webgl_retinomaps analyses
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
sys.argv[5]: recache pycortex
-----------------------------------------------------------------------------------------
Output(s):
Pycortex webgl
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/webgl/
2. run python command
>> python pycortex_webgl_retinomaps.py [main dir] [project] [subject] [group] [recache]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/webgl/
python pycortex_webgl_retinomaps.py /scratch/mszinte/data RetinoMaps sub-01 327 1
python pycortex_webgl_retinomaps.py /scratch/mszinte/data RetinoMaps sub-170k 327 1
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
import cortex

# Personal import
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from pycortex_utils import draw_cortex, set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
recache = sys.argv[5]
if recache == '1': recache = True
else: recache = False

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
analysis_info = settings[0]

prf_task_name = analysis_info['analysis_task_names'][0]
pycortex_subject_template = analysis_info['pycortex_subject_template']

format_ = "170k"

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define directory
pp_dir = "{}/{}/derivatives/pp_data/{}/{}".format(main_dir, project_dir, subject, format_)
rois_datasets_dir = '{}/rois/pycortex/datasets_rois'.format(pp_dir)
css_dataset_dir = "{}/prf/pycortex/datasets_css".format(pp_dir)
intertask_dataset_dir = "{}/intertask/pycortex/datasets_stats".format(pp_dir)

# Define filenames
overlay_fn = "{}/db/{}/overlays_rois-group-mmp.svg".format(cortex_dir, pycortex_subject_template)
rois_datasets_fn = "{}/{}_fmriprep_dct_z-score_rois-group-mmp.hdf".format(rois_datasets_dir, subject )  
css_datasets_fn = "{}/{}_task-{}_fmriprep_dct_z-score_loo-avg_rois-group-mmp_prf-css.hdf".format(css_dataset_dir, subject, prf_task_name)
intertask_datasets_fn = "{}/{}_intertask_Sac-Pur-pRF.hdf".format(intertask_dataset_dir, subject)

# Concatenate filenames
dataset_list_fns = []
dataset_list_fns.append([rois_datasets_fn])
dataset_list_fns.append([css_datasets_fn])
dataset_list_fns.append([intertask_datasets_fn])

# Load datasets and combine them
list_dataset = ''
dataset_num = 0
for dataset_fn_list in dataset_list_fns:
    for dataset_fn in dataset_fn_list:
        dataset_num += 1
        exec("dataset_{} = cortex.load(dataset_fn)".format(dataset_num))
        list_dataset += "dataset_{}=dataset_{}, ".format(dataset_num, dataset_num)
exec("new_dataset = cortex.Dataset({})".format(list_dataset))

views_to_keep = ['ROIs', 'CSS pRF CM', 'CSS pRF eccentricity', 'CSS pRF n', 'CSS pRF polar angle', 'CSS pRF size', 'saccade', 'pursuit', 'vision', 'visuomotor']
new_dataset = cortex.Dataset(**{k: new_dataset.views[k] for k in views_to_keep})

# Make webgl
webgl_dir = "{}/{}/derivatives/webgl/{}/{}".format(main_dir, project_dir, subject, format_)
os.makedirs(webgl_dir, exist_ok=True)
print("Saving: {}".format(webgl_dir))
if os.path.isdir(webgl_dir):os.system("rm -Rfvd {}".format(webgl_dir))    
labels_visible = ('rois')


cortex.webgl.make_static(outpath=webgl_dir,
                         data=new_dataset,
                         overlays_visible=labels_visible,
                         labels_visible=labels_visible, 
                         overlay_file = overlay_fn,
                         title="Project:{}; Subject:{}; Format:{}".format(project_dir, subject, format_),
                         recache=recache, 
                         # curvature_brightness=0.1,  
                         # curvature_contrast=0.8, 
                         # curvature_smoothness=1, 
                         # surface_specularity=0
                         )

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))