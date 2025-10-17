"""
-----------------------------------------------------------------------------------------
prf_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
pRF fit code run by submit_fit.py
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name
sys.argv[2]: input filepath (timeseries)
sys.argv[2]: visual design filepath
sys.argv[3]: model output filepath
sys.argv[4]: timeseries prediction output path
sys.argv[5]: number of processors
-----------------------------------------------------------------------------------------
Output(s):
Nifti image files with fit parameters for a z slice
-----------------------------------------------------------------------------------------
To run :
>> cd to function directory
cd ~/projects/amblyo_prf/analysis_code/postproc/
>> python prf/fit/prf_fit.py [subject] [timeseries] [visual design] 
                     [fit] [prediction] [nb_procs]
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import sys, os
import numpy as np
import datetime
import json
import ipdb
deb = ipdb.set_trace

# MRI analysis imports
from prfpy.rf import *
from prfpy.timecourse import *

#sys.path.insert(0,r'/tmp/prfpy')
import prfpy
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel, Norm_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
import nibabel as nb


# Get inputs
start_time = datetime.datetime.now()

subject = sys.argv[1]

if sys.argv[2].endswith('.nii'):
    input_fn_HCP = sys.argv[2]
elif sys.argv[2].endswith('.gii'):
    input_fn_fsnative = sys.argv[2]

input_vd = sys.argv[3]

if sys.argv[4].endswith('.nii'):
    fit_fn_HCP = sys.argv[4]   
elif sys.argv[4].endswith('.gii'):
    fit_fn_fsnative = sys.argv[4]

if sys.argv[5].endswith('.nii'):
    pred_fn_HCP = sys.argv[5]
elif sys.argv[5].endswith('.gii'):
    pred_fn_fsnative = sys.argv[5]


nb_procs = int(sys.argv[6])

# Analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
screen_size_cm = analysis_info['screen_size_cm']
screen_distance_cm = analysis_info['screen_distance_cm']
TR = analysis_info['TR']
grid_nr = analysis_info['grid_nr']
max_ecc_size = analysis_info['max_ecc_size']

# Get task specific visual design matrix
vdm = np.load(input_vd)

# if input_fn_HCP in locals(): 


# Load HCP 170k data 
data_img_HCP = nb.load(input_fn_HCP)
data_HCP = data_img_HCP.get_fdata()




test_data = data_HCP[1000:1010,:]

data_HCP = test_data

    
fit_mat_gauss_HCP = np.zeros((data_HCP.shape[0],data_HCP.shape[1],6))
pred_mat_gauss_HCP = np.zeros(data_HCP.shape)

fit_mat_DN_HCP = np.zeros((data_HCP.shape[0],data_HCP.shape[1],6))
pred_mat_DN_HCP = np.zeros(data_HCP.shape)

# determine model
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

sizes = max_ecc_size * np.linspace(0.1,1,grid_nr)**2 
eccs = max_ecc_size * np.linspace(0.25,1,grid_nr)**2
polars = np.linspace(0, 2*np.pi, grid_nr)

#create gaussian base model
gauss_model = Iso2DGaussianModel(stimulus=stimulus, filter_predictions=True,filter_type='dc',filter_params={'first_modes_to_remove': 3})


# grid fit
print("Grid fit")
gauss_fitter = Iso2DGaussianFitter(data=test_data, model=gauss_model, n_jobs=nb_procs)
gauss_fitter.grid_fit(ecc_grid=eccs, polar_grid=polars, size_grid=sizes)

   

   

#GAUSSIAN ITERATIVE FIT
print("Iterative fit")
rsq_threshold=0.01
verbose=True
gauss_bounds = [(-1*max_ecc_size, 1*max_ecc_size),  # x
                (-1*max_ecc_size, 1*max_ecc_size),  # y
                (0.1, 1*10),  # prf size
                (0, 1000),  # prf amplitude
                (0, 0)]  # bold baseline
gauss_bounds += [(0,10),(0,0)] #hrf bounds. if want it fixed to some value, specify e.g. (4,4) (0,0)
constraints=None
xtol=1e-4
ftol=1e-4

gauss_fitter.iterative_fit(rsq_threshold=0.01, verbose=verbose, bounds=gauss_bounds,
                 constraints=constraints, 
                 xtol=xtol,
                 ftol=xtol)

fit_fit_gauss_HCP = gauss_fitter.iterative_search_params


#DIVISIVE NORMALISATION MODEL 
print('DN fit')
hrf=[1,0.4,0]

filter_type='dc'
filter_params={"first_modes_to_remove":3,
                        "last_modes_to_remove_percent":0,
                        "window_length":50,
                        "polyorder":3,
                        "highpass":True,
                        "add_mean":True}

filter_predictions=False                                     
normalize_RFs=False

use_previous_gaussian_fitter_hrf=True #if true, will use hrf result from gauss fit at the grid stage instead of doing a grid fit for it

DN_model = Norm_Iso2DGaussianModel(stimulus=stimulus,
                                    hrf=hrf,
                                    filter_predictions=filter_predictions,
                                    filter_type=filter_type,
                                    filter_params=filter_params,                                       
                                    normalize_RFs=normalize_RFs)




DN_fitter = Norm_Iso2DGaussianFitter(data=test_data,
                                model=DN_model,
                                n_jobs=8,
                                previous_gaussian_fitter=gauss_fitter, #need iterative_search params attribute
                                use_previous_gaussian_fitter_hrf=use_previous_gaussian_fitter_hrf)

norm_grid_bounds = [(0,1000),(0,1000)] #only prf amplitudes between 0 and 1000, only neural baseline values between 0 and 1000

DN_fitter.grid_fit(surround_amplitude_grid=np.linspace(0,10,grid_nr),
             surround_size_grid=np.linspace(1,10,grid_nr),
             neural_baseline_grid=np.linspace(0,10,grid_nr),
             surround_baseline_grid=np.linspace(1,10,grid_nr), 
             verbose = False, 
             n_batches = 8, 
             rsq_threshold= 0.01, 
             fixed_grid_baseline=0, 
             grid_bounds=norm_grid_bounds,
             hrf_1_grid=np.linspace(0,10,grid_nr),
             hrf_2_grid=np.linspace(0,0,1),
             ecc_grid=eccs[:grid_nr],
             polar_grid=polars[:grid_nr],
             size_grid=sizes[:grid_nr])


#iterative DN fit
norm_bounds =  [(-1*max_ecc_size, 1*max_ecc_size),  # x
            (-1*max_ecc_size, 1*max_ecc_size),  # y
            (0.1, 1*10),  # prf size
            (0, 1000),  # prf amplitude
            (0, 0),  # bold baseline
            (0, 1000),  # surround amplitude
            (0.1, 3*5),  # surround size
            (0, 1000),  # neural baseline
            (1e-6, 1000)]  # surround baseline

norm_bounds += [(0,10),(0,0)] #hrf bounds
constraints_norm = None

bounds_norm = None

DN_fitter.iterative_fit(rsq_threshold=rsq_threshold, verbose=verbose,
                            bounds=bounds_norm,  #!need to find out why unequal amount of start_param and norm bounds
                            constraints=constraints_norm,
                            xtol=xtol,
                            ftol=ftol)
fit_fit_DN_HCP = DN_fitter.iterative_search_params

# Re-arrange data
for est,vox in enumerate(data_indices_HCP):
    fit_mat_gauss_HCP[vox] = fit_fit_gauss_HCP[est]
    pred_mat_gauss_HCP[vox] = gauss_model.return_prediction(  mu_x=fit_fit_gauss_HCP[est][0], mu_y=fit_fit_gauss_HCP[est][1], size=fit_fit_gauss_HCP[est][2], 
                                                    beta=fit_fit_gauss_HCP[est][3], baseline=fit_fit_gauss_HCP[est][4])

# Re-arrange data
for est,vox in enumerate(data_indices_HCP):
    fit_mat_DN_HCP[vox] = fit_fit_DN_HCP[est]
    pred_mat_DN_HCP[vox] = DN_model.return_prediction(  mu_x=fit_fit_DN_HCP[est][0], mu_y=fit_fit_DN_HCP[est][1], size=fit_fit_DN_HCP[est][2], 
                                                    beta=fit_fit_DN_HCP[est][3], baseline=fit_fit_DN_HCP[est][4])



# export HCP DN
fit_img_DN_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=fit_mat_DN_HCP, header=data_img_HCP.header)
nb.save(fit_img_DN_HCP, fit_fn_HCP)

pred_img_DN_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=pred_mat_DN_HCP, header=data_img_HCP.header)
nb.save(pred_img_DN_HCP, pred_fn_HCP)


# export HCP gauss
fit_img_gauss_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=fit_mat_gauss_HCP, header=data_img_HCP.header)
nb.save(fit_img_gauss_HCP, fit_fn_HCP)

pred_img_gauss_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=pred_mat_gauss_HCP, header=data_img_HCP.header)
nb.save(pred_img_gauss_HCP, pred_fn_HCP)



# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time, end_time=end_time, dur=end_time - start_time))


