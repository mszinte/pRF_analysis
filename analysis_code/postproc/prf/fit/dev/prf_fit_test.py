"""
-----------------------------------------------------------------------------------------
prf_fit_test.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing  grid / iterative no bound and iterative repeated
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: gifti file on which we do the test
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/fit
2. run python command
python prf_fit_test.py [main directory] [project name] [subject name]
-----------------------------------------------------------------------------------------
Exemple:
python prf_fit_test.py /scratch/mszinte/data RetinoMaps sub-02 
                       sub-02_task-pRF_hemi-L_fmriprep_dct_avg_bold.func.gii
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# imports
import sys, os
import numpy as np
import glob
import datetime
import json
from pathlib import Path
import ipdb
deb = ipdb.set_trace

# MRI analysis imports
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter
import nibabel as nb
import cortex

# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
gii_file = sys.argv[4]

# data filenames and folder
pycortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
input_vd = '{}/{}/derivatives/vdm/vdm.npy'.format(main_dir, project_dir)
input_fn_fsnative = '{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_avg/fsnative/{}'.format(
    main_dir, project_dir, subject, gii_file)

prf_fit_test_dir = "{}/{}/derivatives/pp_data/{}/prf/fit/test".format(
    main_dir, project_dir, subject)
os.makedirs(prf_fit_test_dir, exist_ok=True)

# pycortex config file
pycortex_db = "{}/db/".format(pycortex_dir)
pycortex_cm = "{}/colormaps/".format(pycortex_dir)
pycortex_config_file  = cortex.options.usercfg
pycortex_config_file_new = pycortex_config_file[:-4] + '_new.cfg'
pycortex_config_file_old = pycortex_config_file[:-4] + '_old.cfg'

Path(pycortex_config_file_new).touch()
with open(pycortex_config_file, 'r') as fileIn:
    with open(pycortex_config_file_new, 'w') as fileOut:
        for line in fileIn:
            if 'filestore' in line:
                newline = 'filestore=' + pycortex_db
                fileOut.write(newline)
                newline = '\n'
            elif 'colormaps' in line:
                newline = 'colormaps=' + pycortex_cm
                fileOut.write(newline)
                newline = '\n'
            else:
                newline = line
            fileOut.write(newline)
os.rename(pycortex_config_file,pycortex_config_file_old)
os.rename(pycortex_config_file_new, pycortex_config_file)

# analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
screen_size_cm = analysis_info['screen_size_cm']
screen_distance_cm = analysis_info['screen_distance_cm']
TR = analysis_info['TR']
max_ecc_size = analysis_info['max_ecc_size']
gauss_grid_nr = analysis_info['gauss_grid_nr']
css_grid_nr = analysis_info['css_grid_nr']
rsq_iterative_th = analysis_info['rsq_iterative_th']
n_jobs = analysis_info['nb_procs_fit_prf']
n_batches = n_jobs
verbose = False
grid_verbose = False
iterative_verbose = False
rois = ['occ', 'par', 'frt']

# load visual design
vdm = np.load(input_vd)

# determine visual design
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

# Load fsnative data 
data_img_fsnative, data_fsnative = load_surface(input_fn_fsnative)

# Get testing vertices
roi_verts = cortex.get_roi_verts(subject=subject, 
                                 roi=rois, 
                                 mask=True)

vert_mask_occ = roi_verts['occ']
vert_mask_occ = vert_mask_occ[0:data_fsnative.shape[1]]
data_occ = data_fsnative[:, vert_mask_occ].T

vert_mask_par = roi_verts['par']
vert_mask_par = vert_mask_par[0:data_fsnative.shape[1]]
data_par = data_fsnative[:, vert_mask_par].T

vert_mask_frt = roi_verts['frt']
vert_mask_frt = vert_mask_frt[0:data_fsnative.shape[1]]
data_frt = data_fsnative[:, vert_mask_frt].T

data = np.concatenate((data_occ, data_par, data_frt))

# Gaussian model
# --------------

# define gauss model
gauss_model = Iso2DGaussianModel(stimulus=stimulus)

# grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data, 
                                   model=gauss_model, 
                                   n_jobs=n_jobs)

# fit parameters
gauss_params_num = 8
gauss_iterative_bound = True
gauss_iterative_repeat = 10


# gaussian grid fit

# define saving data matrices
gauss_fit_mat_repeat = np.zeros((data.shape[0], gauss_params_num, gauss_iterative_repeat+2))
gauss_pred_mat_repeat = np.zeros((data.shape[0], data.shape[1], gauss_iterative_repeat+2))
gauss_duration_repeat = np.zeros(gauss_iterative_repeat+2)

# defind grid values
ecc_gauss_grid = max_ecc_size * np.linspace(0.25,1.0,gauss_grid_nr)**2
sizes_guass_grid = max_ecc_size * np.linspace(0.1,1.0,gauss_grid_nr)**2
polars_gauss_grid = np.linspace(0.0, 2.0*np.pi, gauss_grid_nr)

# run grid fit
start_time = datetime.datetime.now()
gauss_fitter.grid_fit(ecc_grid=ecc_gauss_grid, 
                      size_grid=sizes_guass_grid, 
                      polar_grid=polars_gauss_grid, 
                      verbose=grid_verbose,
                      n_batches=n_batches)

end_time = datetime.datetime.now()
duration = end_time - start_time
print("Grid fit:\t\t\t{}\nEnd time:\t\t\t{}\nDuration:\t\t\t{}".format(
    start_time, end_time, duration))

# arrange data
gauss_fit = gauss_fitter.gridsearch_params
gauss_fit_mat = np.zeros((data.shape[0],gauss_params_num))
gauss_pred_mat = np.zeros_like(data) 
for est in range(len(data)):
    gauss_fit_mat[est] = gauss_fit[est]
    gauss_pred_mat[est,:] = gauss_model.return_prediction(mu_x=gauss_fit[est][0], 
                                                          mu_y=gauss_fit[est][1], 
                                                          size=gauss_fit[est][2], 
                                                          beta=gauss_fit[est][3], 
                                                          baseline=gauss_fit[est][4],
                                                          hrf_1=gauss_fit[est][5],
                                                          hrf_2=gauss_fit[est][6])

# save across iteration
gauss_fit_mat_repeat[:,:,0] = gauss_fit_mat
gauss_pred_mat_repeat[:,:,0] = gauss_pred_mat
gauss_duration_repeat[0] = duration.total_seconds()

# gaussian iterative fit

# define parameters dict
gauss_iterative_dict = dict(rsq_threshold=rsq_iterative_th, 
                            verbose=iterative_verbose,
                            xtol=1e-4,
                            ftol=1e-4)

# define iterative bounds
gauss_bounds = [(-1.5*max_ecc_size, 1.5*max_ecc_size), # x
                (-1.5*max_ecc_size, 1.5*max_ecc_size), # y
                (0.001, 1.5*max_ecc_size), # prf size
                (0.001, 10.0), # prf amplitude
                (-2.0, 2.0), # bold baseline
                (0.001, 10.0), # hrf_1 bounds
                (0.001, 5.0)] # hrf_2 bounds

# run iterative fit through repetitions
for repeat_num in np.arange(0,gauss_iterative_repeat+1,1):

    # if repeat_num = 0 => no bound
    # if repeat_num >= 1 => bounds
    if repeat_num == 1:
        gauss_iterative_dict['bounds'] = gauss_bounds
    if repeat_num >= 2:
        gauss_iterative_dict['starting_params'] = gauss_fitter.iterative_search_params

    start_time = datetime.datetime.now()
    gauss_fitter.iterative_fit(**gauss_iterative_dict)
    gauss_fit = gauss_fitter.iterative_search_params

    end_time = datetime.datetime.now()
    duration = end_time - start_time
    print("Iterative fit (repeat={}):\t{}\nEnd time:\t\t\t{}\nDuration:\t\t\t{}".format(
        repeat_num, start_time, end_time, duration))
    
    # arrange data
    gauss_fit_mat = np.zeros((data.shape[0],gauss_params_num))
    gauss_pred_mat = np.zeros_like(data) 
    for est in range(len(data)):
        gauss_fit_mat[est] = gauss_fit[est]
        gauss_pred_mat[est,:] = gauss_model.return_prediction(mu_x=gauss_fit[est][0], 
                                                              mu_y=gauss_fit[est][1], 
                                                              size=gauss_fit[est][2], 
                                                              beta=gauss_fit[est][3], 
                                                              baseline=gauss_fit[est][4],
                                                              hrf_1=gauss_fit[est][5],
                                                              hrf_2=gauss_fit[est][6]
                                                             )
        
    # save across iteration
    gauss_fit_mat_repeat[:,:,repeat_num+1] = gauss_fit_mat
    gauss_pred_mat_repeat[:,:,repeat_num+1] = gauss_pred_mat
    gauss_duration_repeat[repeat_num+1] = duration.total_seconds()

# save arrays
np.save(file='{}/{}_gauss_fit_mat_repeat.npy'.format(prf_fit_test_dir, subject),
        arr=gauss_fit_mat_repeat)
np.save(file='{}/{}_gauss_pred_mat_repeat.npy'.format(prf_fit_test_dir, subject),
        arr=gauss_pred_mat_repeat)
np.save(file='{}/{}_gauss_duration_repeat.npy'.format(prf_fit_test_dir, subject),
        arr=gauss_duration_repeat)