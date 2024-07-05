"""
-----------------------------------------------------------------------------------------
pycortex_webgl_css.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create combined webgl per participants for pRF CSS analyses
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
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/webgl/
2. run python command
>> python pycortex_webgl_css.py [main dir] [project] [subject] [group] [recache]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/RetinoMaps/analysis_code/postproc/prf/webgl/
python pycortex_webgl_css.py /scratch/mszinte/data RetinoMaps sub-01 327 1
python pycortex_webgl_css.py /scratch/mszinte/data RetinoMaps sub-170k 327 1
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
import json
import cortex

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
recache = sys.argv[5]
if recache == '1': recache = True
else: recache = False

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
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
    cor_datasets_dir = '{}/corr/pycortex/datasets_inter-run-corr'.format(pp_dir)
    gridfit_datasets_dir = '{}/prf/pycortex/datasets_avg_gauss_gridfit'.format(pp_dir)
    rois_datasets_dir = '{}/rois/pycortex/datasets_rois'.format(pp_dir)
    css_dataset_dir = "{}/prf/pycortex/datasets_loo-avg_css".format(pp_dir)

    # Define filenames
    cor_datasets_fn = []
    cor_datasets_fn.append("{}/{}_task-{}_inter-run-corr.hdf".format(cor_datasets_dir, subject, prf_task_name)) 
    rois_datasets_fn = "{}/{}_task-{}_rois.hdf".format(rois_datasets_dir, subject, prf_task_name)
    gridfit_datasets_fn = "{}/{}_task-{}_avg_gauss_gridfit.hdf".format(gridfit_datasets_dir, subject, prf_task_name)
    css_datasets_fn = "{}/{}_task-{}_loo-avg_css.hdf".format(css_dataset_dir, subject, prf_task_name)

    # Concatenate filenames
    dateset_list_fns = []
    dateset_list_fns.append(cor_datasets_fn)
    dateset_list_fns.append([rois_datasets_fn])
    dateset_list_fns.append([gridfit_datasets_fn])
    dateset_list_fns.append([css_datasets_fn])
    
    # Load datasets and combine them
    list_dataset = ''
    dataset_num = 0
    for dataset_fn_list in dateset_list_fns:
        for dataset_fn in dataset_fn_list:
            dataset_num += 1
            exec("dataset_{} = cortex.load(dataset_fn)".format(dataset_num))
            list_dataset += "dataset_{}=dataset_{}, ".format(dataset_num, dataset_num)
    exec("new_dataset = cortex.Dataset({})".format(list_dataset))
    
    # Make webgl
    webgl_dir = "{}/{}/derivatives/webgl/{}/{}".format(main_dir, project_dir, subject, format_)
    os.makedirs(webgl_dir, exist_ok=True)
    print("Saving: {}".format(webgl_dir))
    if os.path.isdir(webgl_dir):os.system("rm -Rfvd {}".format(webgl_dir))    
    if format_ == 'fsnative': labels_visible = ('rois' )
    else: labels_visible = ('')

    cortex.webgl.make_static(outpath=webgl_dir,
                             data=new_dataset,
                             labels_visible=labels_visible,
                             title="Project:{}; Subject:{}; Format:{}".format(project_dir, subject, format_),
                             recache=recache)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))