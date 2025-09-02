"""
-----------------------------------------------------------------------------------------
pycortex_webgl_intertask.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create combined webgl per participants for intertask analysis 
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
>> python pycortex_webgl_intertask.py [main dir] [project] [subject] [group] [recache]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/webgl/
python pycortex_webgl_intertask.py /scratch/mszinte/data RetinoMaps sub-01 327 1
python pycortex_webgl_intertask.py /scratch/mszinte/data RetinoMaps sub-170k 327 1
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
---------
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
import json
import shutil
import cortex

# Personal import
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
recache = sys.argv[5]
if recache == '1': recache = True
else: recache = False

# Define analysis parameters
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
prf_task_name = analysis_info['prf_task_name']
webapp_dir = analysis_info['webapp_dir']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):

    # Define directory
    pp_dir = "{}/{}/derivatives/pp_data/{}/{}".format(main_dir, project_dir, subject, format_)
    intertask_dataset_dir = '{}/intertask/pycortex/datasets_stats'.format(pp_dir)
    
    # Define filenames
    intertask_datasets_fn = "{}/{}_intertask_Sac_Pur.hdf".format(intertask_dataset_dir, subject)
    
    # Load datasets
    dataset_intertask = cortex.load(intertask_datasets_fn)
    
    # Make webgl
    webgl_dir = "{}/{}/derivatives/webgl/{}/{}".format(main_dir, project_dir, subject, format_)
    os.makedirs(webgl_dir, exist_ok=True)
    print("Saving: {}".format(webgl_dir))
    
    if format_ == 'fsnative': labels_visible = ('rois' )
    else: labels_visible = ('')

    cortex.webgl.make_static(outpath=webgl_dir,
                             data=dataset_intertask,
                             labels_visible=labels_visible,
                             title="Project:{}; Subject:{}; Format:{}".format(project_dir, subject, format_),
                             recache=recache)
    
    # Copy HTML figure file into webgl folder
    html_fig_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/figures'.format(main_dir, project_dir, subject, format_)
    html_fig_fn = '{}/{}_figures_html_Sac_Pur.html'.format(html_fig_dir, subject)
    webgl_html_dir = '{}/figures'.format(webgl_dir)
    os.makedirs(webgl_html_dir, exist_ok=True)
    webgl_html_fn ='{}/{}_figures_html_Sac_Pur.html'.format(webgl_html_dir, subject)
    shutil.copy2(html_fig_fn, webgl_html_fn)
    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
    





