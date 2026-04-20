"""
-----------------------------------------------------------------------------------------
ncsf_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing css fit
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
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python ncsf_fit.py [main directory] [project name] [subject name] 
                     [input file name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
python prf_cssfit.py /scratch/mszinte/data RetinoMaps sub-03 [file path] 32  
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
import sys
import datetime
import numpy as np
import pandas as pd

# MRI analysis imports
import nibabel as nb
from prfpy_csenf.model import CSenFModel
from prfpy_csenf.stimulus import CSenFStimulus
from prfpy_csenf.fit import CSenFFitter
# from prfpy_csenf.rf import * 
# from prfpy_csenf.csenf_plot_functions import *

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings
from pycortex_utils import set_pycortex_config_file
from surface_utils import load_surface, make_surface_image


# Get inputs
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
n_jobs = int(sys.argv[5])

n_batches = n_jobs
verbose = True


# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

TR = analysis_info['TR']
nCSF_ses = 'ses-01'
nCSF_task_name = 'nCSF'
format_ = 'fsnative'
ncsf_grid_nr = 5
rsq_threshold = 0.1 
ncsf_params_num = 10

ncsf_maps_names = ['width_r', 'SFp', 'CSp', 'width_l', 'crf_exp', 'amp_1', 
                   'bold_baseline', 'hrf_1', 'hrf_2', 'r_squared']

n_jobs = 32

sf_minFreq = 0.05
sf_maxFreq = 16
sf_filtNum = 6 

minCont = 0.25
maxCont = 80
contNum = 12

sf_filtCenters = np.concatenate([np.round(np.logspace(np.log10(0.05), np.log10(16), sf_filtNum), 2), [0]])
contValues = np.concatenate([np.logspace(np.log10(minCont), np.log10(maxCont), contNum), [0]])

# Create Stimulus object
# Load events to have stim sequence
events_dir = '{}/{}/{}/{}/func'.format(main_dir, project_dir, subject, nCSF_ses)
events_fn = '{}/{}_{}_task-{}_dir-PA_run-01_events.tsv'.format(events_dir, subject, nCSF_ses, nCSF_task_name)
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

# exlude nan voxel from the analysis 
valid_vertices = ~np.isnan(raw_data).any(axis=0)
valid_vertices_idx = np.where(valid_vertices)[0]
data = raw_data[:,valid_vertices]

# Dermine nCSF model
csenf_model = CSenFModel(stimulus=csenf_stim, 
                         hrf=[1,1,0], 
                         edge_type='CRF')

csenf_fitter = CSenFFitter(data=data.T, 
                           model=csenf_model, 
                           n_jobs=n_jobs)

# nCSF bounds
ncsf_bounds = {
    'width_r' : [0,10],          
    'SFp' : [0, 20],
    'CSp' : [0, 200] ,
    'width_l' : [0.68, 0.68],
    'crf_exp' : [0, 10],
    'amp_1' : [0, 6],
    'bold_baseline' : [-1,1] ,
    'hrf_1' : [1, 1],
    'hrf_2' : [0,0],
}


width_r_grid        = np.linspace(ncsf_bounds['width_r'][0], ncsf_bounds['width_r'][1], ncsf_grid_nr)     
SFp_grid            = np.linspace(ncsf_bounds['SFp'][0], ncsf_bounds['SFp'][1], ncsf_grid_nr)     
CSp_grid            = np.linspace(ncsf_bounds['CSp'][0], ncsf_bounds['CSp'][1], ncsf_grid_nr)
width_l_grid        = np.linspace(ncsf_bounds['width_l'][0], ncsf_bounds['width_l'][1], ncsf_grid_nr)     
crf_exp_grid        = np.linspace(ncsf_bounds['crf_exp'][0], ncsf_bounds['crf_exp'][1], ncsf_grid_nr)
hrf_1_grid = None
hrf_2_grid = None

grid_bounds = [ncsf_bounds['amp_1']]
fixed_grid_baseline=False


bounds_list = [
    (ncsf_bounds['width_r']),     # width_r
    (ncsf_bounds['SFp']),     # SFp
    (ncsf_bounds['CSp']),    # CSp
    (ncsf_bounds['width_l']),     # width_l
    (ncsf_bounds['crf_exp']),     # crf_exp
    (ncsf_bounds['amp_1']),   # amp_1
    (ncsf_bounds['bold_baseline']),      # baseline
    (ncsf_bounds['hrf_1']),      # baseline
    (ncsf_bounds['hrf_2']),      # baseline
]

# Grid fit
csenf_fitter.grid_fit(width_r_grid=width_r_grid, 
                      SFp_grid=SFp_grid, 
                      CSp_grid=CSp_grid, 
                      width_l_grid=width_l_grid, 
                      crf_exp_grid=crf_exp_grid, 
                      hrf_1_grid=hrf_1_grid, 
                      hrf_2_grid=hrf_2_grid, 
                      verbose=True,
                      fixed_grid_baseline=fixed_grid_baseline, 
                      grid_bounds=grid_bounds, 
                      n_batches=n_jobs)

# Iterative fit
csenf_fitter.iterative_fit(
    rsq_threshold = rsq_threshold,
    verbose = False,
    bounds = bounds_list,
    # constraints = csf_constraints,
    # xtol=xtol,   
    # ftol=ftol,           
    )


ncsf_fit = csenf_fitter.iterative_search_params


# rearange result of Gauss model 
#gauss_fit = gauss_fitter.gridsearch_params
ncsf_fit_mat = np.zeros((raw_data.shape[1], ncsf_params_num))
ncsf_pred_mat = np.zeros_like(raw_data) 

for est, vert in enumerate(valid_vertices_idx):
    ncsf_fit_mat[vert] = ncsf_fit[est]
    # ncsf_pred_mat[:,vert] = csenf_model.return_prediction(mu_x=ncsf_fit[est][0], 
    #                                                       mu_y=ncsf_fit[est][1], 
    #                                                       size=ncsf_fit[est][2], 
    #                                                       beta=ncsf_fit[est][3], 
    #                                                       baseline=ncsf_fit[est][4],
    #                                                       hrf_1=ncsf_fit[est][5],
    #                                                       hrf_2=ncsf_fit[est][6])


ncsf_fit_mat = np.where(ncsf_fit_mat == 0, np.nan, ncsf_fit_mat)
ncsf_pred_mat = np.where(ncsf_pred_mat == 0, np.nan, ncsf_pred_mat)

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
img_ncsf_fit_mat = make_surface_image(data=ncsf_fit_mat.T, source_img=img, maps_names=ncsf_maps_names)
nb.save(img_ncsf_fit_mat,'{}/{}'.format(ncsf_fit_dir, ncsf_fit_fn)) 

# export pred
img_ncsf_pred_mat = make_surface_image(data=ncsf_pred_mat, source_img=img)
nb.save(img_ncsf_pred_mat,'{}/{}'.format(ncsf_fit_dir, ncsf_pred_fn)) 


























