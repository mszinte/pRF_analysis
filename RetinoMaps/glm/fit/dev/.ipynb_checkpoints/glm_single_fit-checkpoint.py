"""
-----------------------------------------------------------------------------------------
glm_single_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make GLM fit with GLM signle
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# GLM  prediction
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/glm/fit/dev
2. run python command
python glm_fit.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/glm/fit/dev
python glm_single_fit.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
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
import re
import sys
import json
import glob
import datetime
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt

# GLMsingle
from glmsingle.glmsingle import GLM_single
from glmsingle.glmsingle import getcanonicalhrflibrary
from glmsingle.utils.squish import squish

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from glm_utils import make_design, eventsMatrix
from surface_utils import make_surface_image
from pycortex_utils import data_from_rois, set_pycortex_config_file

# Start counting the elapsed time for code execution
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
TR = analysis_info['TR']
TRs = analysis_info['TRs']
tasks = analysis_info['task_glm']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
func_session = analysis_info['func_session'][0]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

for format_, extension in zip(formats, extensions):
    print(format_)
    # make folders
    glm_dir = '{}/{}/derivatives/pp_data/{}/{}/glm/glm_single_fit'.format(
        main_dir, project_dir, subject, format_)
    os.makedirs(glm_dir, exist_ok=True)
    
    if format_ == 'fsnative': rois = analysis_info['rois']
    elif format_ == '170k': rois = analysis_info['mmp_rois']

    for task in tasks:
        print(task)
        # Sub-01 exeption 
        if subject == 'sub-01' and task in ('SacLoc', 'PurLoc'): func_session = 'ses-01'
        elif subject == 'sub-01' and task in ('SacVELoc', 'PurVELoc'): func_session = 'ses-03'

        # prepoc files name
        preproc_fns = glob.glob('{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct_loo_avg/*task-{}*dct_avg_loo*.{}'.format(
            main_dir, project_dir, subject, format_, task, extension))
        
        # Split filtered files  depending of their nature
        preproc_fsnative_hemi_L, preproc_fsnative_hemi_R = [], []
        for subtype in preproc_fns:
            if "hemi-L" in subtype:
                preproc_fsnative_hemi_L.append(subtype)
            elif "hemi-R" in subtype:
                preproc_fsnative_hemi_R.append(subtype)

                
        preproc_files_list = [preproc_fsnative_hemi_L,
                              preproc_fsnative_hemi_R]

        
        for preproc_files in preproc_files_list:
            data_multi_run = []
            design_multi_run = []
            for preproc_fn in preproc_files :
                # match = re.search(r'_loo-(\d+)_', preproc_fn)
                # loo_num = 'loo-{}'.format(match.group(1))
            
                # find the events and confounds files 
                event_dir = '{}/{}/{}/{}/func/'.format(
                    main_dir, project_dir, subject, func_session)
                
                con_dir = '{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func'.format(
                    main_dir, project_dir, subject, func_session)
        
                # Find the event files 
                event_file = glob.glob("{}/{}_{}_task-{}_run-01_events.tsv".format(
                    event_dir, subject, func_session, task))
    
                # make the designe matrixe  
                events = eventsMatrix(design_file=event_file[0], task=task, tr=TR)
                
                frame_times = np.arange(TRs) * TR
                design_matrix = make_design(events=events, tr=TR, n_times=TRs)
                design_matrix_reshaped = design_matrix[:,1].reshape(-1,1)
                
            
                
                
                # #  Save the designe matrix 
                # dm_dir = '{}/{}/derivatives/pp_data/{}/{}/glm/designe_matrix'.format(
                #     main_dir, project_dir, subject, format_)
                # os.makedirs(dm_dir, exist_ok=True)
                # design_matrix.to_csv('{}/{}_task-{}_designe_matrix.tsv'.format(
                #     dm_dir, subject, task), sep="\t", na_rep='NaN', index=False)
                
                # print('Save {}/{}_task-{}_designe_matrix.pdf'.format(
                #     dm_dir, subject, task))
                # plt.figure()
                # plot_design_matrix(design_matrix)
                # plt.savefig('{}/{}_task-{}_designe_matrix.pdf'.format(
                #     dm_dir, subject, task))
                # plt.close()
    
                # make glm output filenames
                glm_pred_fn = preproc_fn.split('/')[-1].replace('bold', 'glm-pred') 
                glm_fit_fn = preproc_fn.split('/')[-1].replace('bold', 'glm-fit') 
    
                # Load data
                preproc_img, preproc_data_brain, preproc_data_rois, roi_idx = data_from_rois(
                    fn=preproc_fn,  subject=subject, rois=rois)
                
                # Make the data and design multi run list 
                data_multi_run.append(preproc_data_rois.T)
                design_multi_run.append(design_matrix_reshaped)
            
            
            # Initialize the model
            opt = dict()
            opt['wantlibrary'] = 1
            opt['wantglmdenoise'] = 1
            opt['wantfracridge'] = 1
            opt['wantfileoutputs'] = [1,1,1,1]
            opt['wantmemoryoutputs'] = [1,1,1,1]
            glmsingle_obj = GLM_single(opt)
            
            
            
            
            # Run the glm
            outputdir = '{}/output'.format(glm_dir)
            os.makedirs(outputdir, exist_ok=True)
            figdir = '{}/figure'.format(glm_dir)
            os.makedirs(figdir, exist_ok=True)
            
            stimdur= 38.4
            
            results_glmsingle = glmsingle_obj.fit(design=design_multi_run, 
                                                  data=data_multi_run, 
                                                  stimdur=stimdur, 
                                                  tr=TR, 
                                                  outputdir=outputdir, 
                                                  figuredir=figdir)
            
            deb()
            # make the model prediction 
            # Import the designSINGLE
            designSINGLE = np.load('{}/DESIGNINFO.npy'.format(outputdir), allow_pickle=True).item()['designSINGLE']
            
            # acces to the library of HRF 
            hrflib = getcanonicalhrflibrary(stimdur,TR).transpose()
            
            
            for n_run in range(len(data_multi_run)):
                designSINGLE_run = designSINGLE[n_run]
                betatemp_all = np.zeros_like(betas_all)
                glmsingle_pred = np.zeros(data_multi_run[n_run].shape)
                
                
            
            
            
            
            
            
            
            
            
            
            
            for vert, roi_idx in enumerate(roi_idx):
                glm_pred_brain[:,roi_idx] = glm_pred_rois[:,vert]
        
            # export pred
            print('Save {}/{}'.format(glm_dir, glm_pred_fn))
            pred_img = make_surface_image(data=glm_pred_brain, 
                                          source_img=preproc_img)
            nb.save(pred_img,'{}/{}'.format(glm_dir, glm_pred_fn)) 
            

                    
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))
       
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))
