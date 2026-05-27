"""
-----------------------------------------------------------------------------------------
ncsf_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Fit nCSF model to data
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: number of jobs 
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
2. run python command
python ncsf_fit.py [main directory] [project name] [subject name] 
                     [input file name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
python ncsf_fit.py /scratch/mszinte/data nCSF sub-01 /scratch/mszinte/data/nCSF/derivatives/pp_data/sub-01/fsnative/func/fmriprep_dct_z-score_avg/sub-01_task-nCSF_hemi-R_fmriprep_dct_z-score_avg_bold.func.gii 16  
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
import datetime
import numpy as np
import pandas as pd

# MRI analysis imports
import nibabel as nb
from prfpy_csenf.fit import CSenFFitter
from prfpy_csenf.model import CSenFModel
from prfpy_csenf.stimulus import CSenFStimulus
# from prfpy_csenf.rf import * 
from prfpy_csenf.csenf_plot_functions import ncsf_calculate_sfmax, ncsf_calculate_aulcsf

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from maths_utils import r2_score_surf
from settings_utils import load_settings
from pycortex_utils import set_pycortex_config_file, get_rois
from surface_utils import load_surface, make_surface_image

start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
n_jobs = int(sys.argv[5])

n_batches = n_jobs
verbose = False

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
nCSF_settings_path = os.path.join(base_dir, project_dir, "nCSF-analysis.yml")
settings = load_settings([settings_path, prf_settings_path, nCSF_settings_path])
settings = load_settings([settings_path, nCSF_settings_path, prf_settings_path])
analysis_info = settings[0]

TR = analysis_info['TR']
nCSF_session = analysis_info['nCSF_session']
rois_methods = analysis_info['rois_methods']
ncsf_grid_nr = analysis_info['ncsf_grid_nr']
nCSF_task_name = analysis_info['nCSF_task_name']
run_ncsf_fit_on_rois = analysis_info['run_ncsf_fit_on_rois']
rsq_threshold = analysis_info['ncsf_iterative_rsq_threshold']
pycortex_subject_template = analysis_info['pycortex_subject_template']

sf_minFreq = analysis_info['sf_minFreq']
sf_maxFreq = analysis_info['sf_maxFreq']
sf_filtNum = analysis_info['sf_filtNum']
minCont = analysis_info['minCont']
maxCont = analysis_info['maxCont']
contNum = analysis_info['contNum']

sf_filtCenters = np.concatenate([np.round(np.logspace(np.log10(0.05), np.log10(16), sf_filtNum), 2), [0]])
contValues = np.concatenate([np.logspace(np.log10(minCont), np.log10(maxCont), contNum), [0]])

maps_names_ncsf = analysis_info['maps_names_ncsf']
for idx, col_name in enumerate(maps_names_ncsf):
    exec("{}_idx = idx".format(col_name))

# Create Stimulus object
# Load events to have stim sequence
events_dir = '{}/{}/{}/{}/func'.format(main_dir, project_dir, subject, nCSF_session)
events_fn = '{}/{}_{}_task-{}_dir-PA_run-01_events.tsv'.format(events_dir, subject, nCSF_session, nCSF_task_name)
events_df = pd.read_table(events_fn, sep="\t")

# Create mapping betwen SF, MC numbers and real values
mapping_SF = {i+1: sf_filtCenters[i] for i in range(len(sf_filtCenters))}
mapping_MC = {i+1: contValues[i] for i in range(len(contValues))}

events_df['spatial_frequency'] = events_df['spatial_frequency'].replace(mapping_SF)
events_df['michelson_contrast'] = events_df['michelson_contrast'].replace(mapping_MC)

sfs_seq = np.array(events_df['spatial_frequency'])
con_seq = np.array(events_df['michelson_contrast'])

csenf_stim = CSenFStimulus(
    SF_seq  = sfs_seq, # np array, 1D 
    CON_seq = con_seq, # np array, 1D 
    TR      = TR,
    discrete_levels=True, # To Check !!!!!! 
    )

# load data
img, raw_data = load_surface(fn=input_fn)

# Extract data from rois if partial fitting
if run_ncsf_fit_on_rois : 
    if input_fn.endswith('.nii'):
       format_ = '170k'
       pycortex_subject =  pycortex_subject_template
       hemi = None
       
    elif input_fn.endswith('.gii'):
       format_ = 'fsnative'
       pycortex_subject =  subject
       match = re.search(r'hemi-[LR]', input_fn)
       hemi = match.group()
    
    rois_method_format = rois_methods[format_][0]
    if rois_method_format == 'rois-drawn':
        rois = analysis_info[rois_method_format]
    elif rois_method_format == 'rois-group-mmp':
        rois = list(analysis_info[rois_method_format].keys())
   
    print('Running fit only on: {}'.format(rois))
    roi_verts_dict = get_rois(pycortex_subject, 
                              surf_format=format_, 
                              rois_type=rois_method_format,
                              mask=False, 
                              rois=rois, 
                              hemis=hemi)
    
    valid_vertices = ~np.isnan(raw_data).any(axis=0)
    roi_vertices_idx = np.unique(np.concatenate(list(roi_verts_dict.values())))
    valid_vertices_idx = roi_vertices_idx[valid_vertices[roi_vertices_idx]]
    data = raw_data[:,valid_vertices_idx]

else:
    len_data = raw_data.shape[1]
    valid_mask = ~np.isnan(raw_data).any(axis=0)
    valid_vertices_idx = np.where(valid_mask)[0]
    data = raw_data[:, valid_mask]


# Dermine nCSF model
csenf_model = CSenFModel(stimulus=csenf_stim, 
                         hrf=[1,1,0], 
                         edge_type='CRF')

csenf_fitter = CSenFFitter(data=data.T, 
                           model=csenf_model, 
                           n_jobs=n_jobs)

ncsf_bounds = {
    'width_r'       : [0, 1.5],          
    'SFp'           : [sf_filtCenters[1], sf_maxFreq],
    'CSp'           : [100/maxCont, 200],
    'width_l'       : [0.5, 1],
    'crf_exp'       : [0, 4],
    'amp_1'         : [0, 6],
    'bold_baseline' : [-1, 1],
    'hrf_1'         : [1, 1],
    'hrf_2'         : [0, 0],
}

width_r_grid  = np.linspace(max(ncsf_bounds['width_r'][0], 0.1), ncsf_bounds['width_r'][1], ncsf_grid_nr)
SFp_grid      = np.logspace(np.log10(ncsf_bounds['SFp'][0]), np.log10(ncsf_bounds['SFp'][1]), ncsf_grid_nr)
CSp_grid      = np.logspace(np.log10(ncsf_bounds['CSp'][0]), np.log10(ncsf_bounds['CSp'][1]), ncsf_grid_nr)
width_l_grid  = np.linspace(ncsf_bounds['width_l'][0], ncsf_bounds['width_l'][1], ncsf_grid_nr)
crf_exp_grid  = np.linspace(max(ncsf_bounds['crf_exp'][0], 0.1), ncsf_bounds['crf_exp'][1], ncsf_grid_nr)
hrf_1_grid    = np.zeros(len(width_r_grid))
hrf_2_grid    = np.zeros(len(width_r_grid))

grid_bounds = [ncsf_bounds['amp_1']]
fixed_grid_baseline=False

bounds_list = [
    (ncsf_bounds['width_r']),
    (ncsf_bounds['SFp']),
    (ncsf_bounds['CSp']),
    (ncsf_bounds['width_l']),
    (ncsf_bounds['crf_exp']),
    (ncsf_bounds['amp_1']),
    (ncsf_bounds['bold_baseline']),
    (ncsf_bounds['hrf_1']),
    (ncsf_bounds['hrf_2']),
]

# Grid fit
csenf_fitter.grid_fit(width_r_grid=width_r_grid, 
                      SFp_grid=SFp_grid, 
                      CSp_grid=CSp_grid, 
                      width_l_grid=width_l_grid, 
                      crf_exp_grid=crf_exp_grid, 
                      hrf_1_grid=hrf_1_grid, 
                      hrf_2_grid=hrf_2_grid, 
                      verbose=verbose,
                      fixed_grid_baseline=fixed_grid_baseline, 
                      grid_bounds=grid_bounds, 
                      n_batches=n_jobs)

# Iterative fit
csenf_fitter.iterative_fit(
    rsq_threshold = rsq_threshold,
    verbose = verbose,
    bounds = bounds_list,
    # constraints = csf_constraints,
    # xtol=xtol,   
    # ftol=ftol,           
    )

ncsf_fit = csenf_fitter.iterative_search_params

# compute AUC and normalise auc and add them to ncsf_fit
auc = ncsf_calculate_aulcsf(width_r=ncsf_fit[:, width_r_idx], 
                            SFp=ncsf_fit[:, SFp_idx], 
                            CSp=ncsf_fit[:, CSp_idx], 
                            width_l=ncsf_fit[:, width_l_idx], 
                            normalize_AUC=False, 
                            SF_levels=sf_filtCenters[:-1]) 

normalise_auc = ncsf_calculate_aulcsf(width_r=ncsf_fit[:, width_r_idx], 
                                      SFp=ncsf_fit[:, SFp_idx], 
                                      CSp=ncsf_fit[:, CSp_idx], 
                                      width_l=ncsf_fit[:, width_l_idx], 
                                      normalize_AUC=True, 
                                      SF_levels=sf_filtCenters[:-1]) 

SFmax = ncsf_calculate_sfmax(width_r=ncsf_fit[:, width_r_idx], 
                             SFp=ncsf_fit[:, SFp_idx], 
                             CSp=ncsf_fit[:, CSp_idx], 
                             max_sfmax=50)

ncsf_fit_extend = np.hstack([ncsf_fit, auc.reshape(-1,1), 
                             normalise_auc.reshape(-1,1), 
                             SFmax.reshape(-1,1)])

# rearange result of nCSF model 
#gauss_fit = gauss_fitter.gridsearch_params
ncsf_fit_mat = np.zeros((raw_data.shape[1], ncsf_fit_extend.shape[1]))
ncsf_pred_mat = np.zeros_like(raw_data) 

for est, vert in enumerate(valid_vertices_idx):
    ncsf_fit_mat[vert] = ncsf_fit_extend[est]
    ncsf_pred_mat[:,vert] = csenf_model.return_prediction(width_r=ncsf_fit[est][width_r_idx], 
                                                          SFp=ncsf_fit[est][SFp_idx], 
                                                          CSp=ncsf_fit[est][CSp_idx], 
                                                          width_l=ncsf_fit[est][width_l_idx], 
                                                          crf_exp=ncsf_fit[est][crf_exp_idx],
                                                          beta=ncsf_fit[est][amp_1_idx],
                                                          baseline=ncsf_fit[est][bold_baseline_idx],
                                                          hrf_1=ncsf_fit[est][hrf_1_idx],
                                                          hrf_2=ncsf_fit[est][hrf_2_idx])

# ncsf_fit_mat = np.where(ncsf_fit_mat == 0, np.nan, ncsf_fit_mat)
# ncsf_pred_mat = np.where(ncsf_pred_mat == 0, np.nan, ncsf_pred_mat)

# Compute LOO r2 
if 'loo-avg' in input_fn:
    # Compute loo r2
    loo_bold_fn = input_fn.replace('loo-avg-', 'loo-')
    loo_img, loo_bold = load_surface(fn=loo_bold_fn)
    loo_r2 = r2_score_surf(bold_signal=loo_bold, model_prediction=ncsf_pred_mat)
    
    # Add loo r2 css_fit_mat
    ncsf_fit_mat = np.column_stack((ncsf_fit_mat, loo_r2))

    # Export data
    maps_names_ncsf = maps_names_ncsf + ['ncsf_loo_rsq']

#export data from gauss model fit

# Define directories
if input_fn.endswith('.nii'):
    ncsf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/ncsf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(ncsf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    ncsf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/ncsf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(ncsf_fit_dir, exist_ok=True)

ncsf_fit_fn  = input_fn.split('/')[-1]
ncsf_fit_fn = ncsf_fit_fn.replace('bold', 'ncsf_fit')

ncsf_pred_fn = input_fn.split('/')[-1]
ncsf_pred_fn = ncsf_pred_fn.replace('bold', 'ncsf_pred')

# export fit
img_ncsf_fit_mat = make_surface_image(data=ncsf_fit_mat.T, source_img=img, maps_names=maps_names_ncsf)
nb.save(img_ncsf_fit_mat,'{}/{}'.format(ncsf_fit_dir, ncsf_fit_fn)) 

# export pred
img_ncsf_pred_mat = make_surface_image(data=ncsf_pred_mat, source_img=img)
nb.save(img_ncsf_pred_mat,'{}/{}'.format(ncsf_fit_dir, ncsf_pred_fn)) 

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
start_time=start_time, end_time=end_time, dur=end_time - start_time))
























