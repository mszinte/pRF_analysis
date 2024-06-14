def fit2deriv(fit_array, model,is_loo_r2=False):
    """
    Compute pRF derivatives out of fitting output and predictions

    Parameters
    ----------
    fit_array: fit parameters 2D array (fit parameter, vertex)
    model: model use for the fit ('gauss','dn','css')
    
    Returns
    -------
    deriv_array: 2D array with pRF derivatives

    stucture output:
    columns: 1->size of input
    rows: derivatives parameters
    
    """

    # Imports
    # -------
    # General imports
    import numpy as np


    # Compute derived measures from prfs/pmfs
    # ---------------------------------------
    # get data index
    
    if model == 'gauss':
        x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, \
        hrf_1_idx, hrf_2_idx, rsq_idx  = 0, 1, 2, 3, 4, 5, 6, 7
        n_params = 8
    
    elif model == 'dn':
        x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, srf_amplitude_idx, \
        srf_size_idx, neural_baseline_idx, surround_baseline_idx, hrf_1_idx, \
        hrf_2_idx,rsq_idx  = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 
        
        # srf_amplitude
        srf_amplitude = fit_array[srf_amplitude_idx,...]
        
        # srf_size
        srf_size = fit_array[srf_size_idx,...]
        
        # neural_baseline 
        neural_baseline = fit_array[neural_baseline_idx,...]
        
        # surround_baseline
        surround_baseline = fit_array[surround_baseline_idx,...]
        
        n_params = 12
    
    elif model == 'css':
        x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, n_idx, hrf_1_idx, hrf_2_idx,rsq_idx \
        = 0, 1, 2, 3, 4, 5, 6, 7, 8
        n = fit_array[n_idx,...]
        
        # exponentiel param 
        n_params = 9
    
    if is_loo_r2 :
        loo_r2_idx = -1
        n_params += 1
        loo_r2 = fit_array[loo_r2_idx,...]
            

    # change to nan empty voxels
    fit_array[fit_array[...,rsq_idx]==0,...] = np.nan
    
    # r-square
    rsq = fit_array[rsq_idx,...]

    # eccentricity
    ecc = np.nan_to_num(np.sqrt(fit_array[x_idx,...]**2 + fit_array[y_idx,...]**2))

    # polar angle
    complex_polar = fit_array[x_idx,...] + 1j * fit_array[y_idx,...]
    normed_polar = complex_polar / np.abs(complex_polar)
    polar_real = np.real(normed_polar)
    polar_imag = np.imag(normed_polar)
    
    # size
    size_ = fit_array[sigma_idx,...].astype(np.float64)
    size_[size_<1e-4] = 1e-4

    # amplitude
    amp = fit_array[beta_idx,...]
    
    # baseline
    baseline = fit_array[baseline_idx,...]

    # x
    x = fit_array[x_idx,...]

    # y
    y = fit_array[y_idx,...]
    
    # hrf_1
    hrf_1 = fit_array[hrf_1_idx,...]
    
    # hrf_2
    hrf_2 = fit_array[hrf_2_idx,...]
    
    
    # Save results

    deriv_array = np.zeros((n_params + 3,fit_array.shape[1],))*np.nan

    deriv_array[0,...] = rsq
    deriv_array[1,...] = ecc
    deriv_array[2,...] = polar_real
    deriv_array[3,...] = polar_imag
    deriv_array[4,...] = size_
    deriv_array[5,...] = amp
    deriv_array[6,...] = baseline
    deriv_array[7,...] = x
    deriv_array[8,...] = y
    deriv_array[9,...] = hrf_1
    deriv_array[10,...] = hrf_2
    
    if model == 'dn':
        deriv_array[11,...] = srf_amplitude
        deriv_array[12,...] = srf_size
        deriv_array[13,...] = neural_baseline
        deriv_array[14,...] = surround_baseline
    
    if model == 'css':
        deriv_array[11,...] = n
        
    # Include leave-one-out R2 if requested
    if is_loo_r2:
        deriv_array[-1, ...] = loo_r2

    
    

    deriv_array = deriv_array.astype(np.float32)

    return deriv_array










