# Debug
import ipdb
deb = ipdb.set_trace
def weighted_regression(x_reg, y_reg, weight_reg, model):
    """
    Function to compute regression parameter weighted by a matrix (e.g. r2 value),
    where the regression model is y = 1/(cx) + d.

    Parameters
    ----------
    x_reg : array (1D)
        x values to regress
    y_reg : array
        y values to regress
    weight_reg : array (1D) 
        weight values (0 to 1) for weighted regression
    model : str
        Type of regression model, either 'pcm' for the original model or 'linear' for a linear model.

    Returns
    -------
    coef_reg : float or array
        regression coefficient(s)
    intercept_reg : float or str
        regression intercept or a string indicating no intercept (for linear model)
    """

    import numpy as np
    from scipy.optimize import curve_fit
    from sklearn import linear_model
    import ipdb
    deb = ipdb.set_trace
    
    x_reg = np.array(x_reg)
    y_reg = np.array(y_reg)
    
    weight_reg = np.array(weight_reg)

    # Filter out NaN values
    x_reg_nan = x_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
    y_reg_nan = y_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
    weight_reg_nan = weight_reg[~np.isnan(weight_reg)]
    
    if model == 'pcm':
        # Define the model function
        def model_function(x, c, d):
            return 1 / (c * x + d)

        if weight_reg_nan.size >= 2:
            # Perform curve fitting
            params, _ = curve_fit(model_function, x_reg_nan, y_reg_nan, sigma=weight_reg_nan, maxfev=10000)
            c, d = params
        else:
            c, d = np.nan, np.nan
        return c, d

    elif model == 'linear':
        if weight_reg_nan.size >= 2:
            regr = linear_model.LinearRegression()
            
            # Filter out NaN values
            x_reg_nan = x_reg_nan.reshape(-1, 1)
            y_reg_nan = y_reg_nan.reshape(-1, 1)
            
            regr.fit(x_reg_nan, y_reg_nan, sample_weight=weight_reg_nan)
            coef_reg, intercept_reg = regr.coef_[0][0], regr.intercept_[0]
        else: 
            coef_reg, intercept_reg = np.nan, np.nan
        return coef_reg, intercept_reg
    else:
        raise ValueError("Invalid model type. Supported models are 'pcm' and 'linear'.")

def weighted_nan_mean(data, weights):
    """
    Calculate the weighted mean of an array, ignoring NaN values.

    Parameters:
    data (np.ndarray): Array of data points, may contain NaN values.
    weights (np.ndarray): Array of weights corresponding to the data points.

    Returns:
    float: The weighted mean of the data points, ignoring NaN values.
    """
    import numpy as np 
    # Mask NaN values in the data
    mask = ~np.isnan(data)
    
    # Apply the mask to data and weights
    masked_data = data[mask]
    masked_weights = weights[mask]
    
    # Calculate the weighted mean
    mean = np.sum(masked_data * masked_weights) / np.sum(masked_weights)
    return mean

def weighted_nan_median(data, weights):
    """
    Calculate the weighted median of a data array, ignoring NaN values.

    Parameters:
    data (pd.Series, pd.DataFrame, np.ndarray): Data points, may contain NaN values.
    weights (pd.Series, pd.DataFrame, np.ndarray): Weights corresponding to the data points.

    Returns:
    float: The weighted median of the data points, ignoring NaN values.
           Returns NaN if the cumulative weights are not defined.
    """
    import numpy as np
    import pandas as pd

    # Convert data and weights to pandas Series if they are numpy arrays
    if isinstance(data, np.ndarray):
        data = pd.Series(data)
    if isinstance(weights, np.ndarray):
        weights = pd.Series(weights)
        
    # If data and weights are DataFrames, ensure they have a single column
    if isinstance(data, pd.DataFrame):
        if data.shape[1] != 1:
            raise ValueError("DataFrame data must have exactly one column")
        data = data.iloc[:, 0]
        
    if isinstance(weights, pd.DataFrame):
        if weights.shape[1] != 1:
            raise ValueError("DataFrame weights must have exactly one column")
        weights = weights.iloc[:, 0]

    # Mask NaN values in the data
    mask = ~data.isna()

    # Apply the mask to data and weights
    masked_data = data[mask].reset_index(drop=True)
    masked_weights = weights[mask].reset_index(drop=True)

    # Check if there are no valid data points
    if masked_data.size == 0 or masked_weights.size == 0:
        return np.nan

    # Sort the data and corresponding weights
    masked_data = pd.to_numeric(masked_data, errors='coerce').dropna()
    sorted_indices = np.argsort(masked_data)
    sorted_data = masked_data.iloc[sorted_indices].reset_index(drop=True)
    sorted_weights = masked_weights.iloc[sorted_indices].reset_index(drop=True)

    # Calculate the cumulative sum of weights
    cumulative_weights = np.cumsum(sorted_weights)

    # Check if cumulative_weights is defined and has elements
    if cumulative_weights.size == 0:
        return np.nan

    # Find the median position
    median_weight = cumulative_weights.iloc[-1] / 2.0

    # Find the index where the cumulative weight crosses the median weight
    median_index = np.searchsorted(cumulative_weights, median_weight)

    return sorted_data.iloc[median_index]

def weighted_nan_percentile(data, weights, percentile):
    """
    Calculate the weighted percentile of an array or a pandas Series, ignoring NaN values.

    Parameters:
    data (np.ndarray or pd.Series): Array or pandas Series of data points, may contain NaN values.
    weights (np.ndarray or pd.Series): Array or pandas Series of weights corresponding to the data points.
    percentile (float): Percentile to compute, between 0 and 100.

    Returns:
    float: The weighted percentile of the data points, ignoring NaN values.
           Returns NaN if the cumulative weights are not defined.
    """
    import numpy as np
    import pandas as pd
    
    # Convert pandas Series to numpy array if needed
    if isinstance(data, pd.Series):
        data = data.values
    if isinstance(weights, pd.Series):
        weights = weights.values
    
    # Mask NaN values in the data
    mask = ~np.isnan(data)
    
    # Apply the mask to data and weights
    masked_data = data[mask]
    masked_weights = weights[mask]
    
    # Check if there are no valid data points
    if masked_data.size == 0 or masked_weights.size == 0:
        return np.nan
    
    # Sort the data and corresponding weights
    sorted_indices = np.argsort(masked_data)
    sorted_data = masked_data[sorted_indices]
    sorted_weights = masked_weights[sorted_indices]
    
    # Calculate the cumulative sum of weights
    cumulative_weights = np.cumsum(sorted_weights)
    
    # Check if cumulative_weights is defined and has elements
    if cumulative_weights.size == 0:
        return np.nan
    
    # Calculate the percentile weight
    percentile_weight = percentile / 100.0 * cumulative_weights[-1]
    
    # Find the index where the cumulative weight crosses the percentile weight
    percentile_index = np.searchsorted(cumulative_weights, percentile_weight)
    
    return float(sorted_data[percentile_index])
        
def gaus_2d(gauss_x, gauss_y, gauss_sd, screen_side, grain=200):
    """
    Generate 2D gaussian mesh
    
    Parameters
    ----------
    gauss_x : mean x gaussian parameter in dva (e.g. 1 dva)
    gauss_y : mean y gaussian parameter in dva (e.g. 1 dva)
    gauss_sd : sd gaussian parameter in dva (e.g. 1 dva)
    screen_side : mesh screen side (square) im dva (e.g. 20 dva from -10 to 10 dva)
    grain : grain resolution of the mesh in pixels (default = 100 pixels)
    
    Returns
    -------
    x : linspace x of the mesh
    y : linspace x of the mesh
    z : mesh_z values (to plot)
    
    """
    import numpy as np
    x = np.linspace(-screen_side/2, screen_side/2, grain)
    y = np.linspace(-screen_side/2, screen_side/2, grain)
    mesh_x, mesh_y = np.meshgrid(x,y) 
    
    gauss_z = 1./(2.*np.pi*gauss_sd*gauss_sd)*np.exp(-((mesh_x-gauss_x)**2./(2.*gauss_sd**2.)+(mesh_y-gauss_y)**2./(2.*gauss_sd**2.)))
    return x, y, gauss_z

def bootstrap_ci_median(data, n_bootstrap=1000, ci_level=0.95):
    import numpy as np
    n = len(data)
    bootstrap_samples = np.random.choice(data, size=(n_bootstrap, n), replace=True)
    medians = np.nanmedian(bootstrap_samples, axis=1)
    lower_ci = np.percentile(medians, (1 - ci_level) / 2 * 100)
    upper_ci = np.percentile(medians, (1 + ci_level) / 2 * 100)
    return lower_ci, upper_ci

def bootstrap_ci_mean(data, n_bootstrap=1000, ci_level=0.95):
    import numpy as np
    n = len(data)
    bootstrap_samples = np.random.choice(data, size=(n_bootstrap, n), replace=True)
    means = np.mean(bootstrap_samples, axis=1)
    lower_ci = np.percentile(means, (1 - ci_level) / 2 * 100)
    upper_ci = np.percentile(means, (1 + ci_level) / 2 * 100)
    return lower_ci, upper_ci

def r2_score_surf(bold_signal, model_prediction):
    """
    Compute r2 between bold signal and model. The gestion of nan values 
    is down with created a non nan mask on the model prediction 

    Parameters
    ----------
    bold_signal: bold signal in 2-dimensional np.array (time, vertex)
    model_prediction: model prediction in 2-dimensional np.array (time, vertex)
    
    Returns
    -------
    r2_scores: the R2 score for each vertex
    """
    import numpy as np
    from sklearn.metrics import r2_score
    
    # Check for NaN values in both bold_signal and model_prediction
    nan_mask = np.isnan(model_prediction).any(axis=0) | np.isnan(bold_signal).any(axis=0)
    valid_vertices = ~nan_mask
    
    # Set R2 scores for vertices with NaN values to NaN
    r2_scores = np.full_like(nan_mask, np.nan, dtype=float)
    
    # Compute R2 scores for vertices without NaN values
    r2_scores[valid_vertices] = r2_score(bold_signal[:, valid_vertices], model_prediction[:, valid_vertices], multioutput='raw_values')
    
    return r2_scores

def linear_regression_surf(bold_signal, model_prediction, correction=None, alpha=None, use_fisher=False):
    """
    Perform linear regression analysis between model predictions and BOLD signals across vertices.

    Parameters:
    bold_signal (numpy.ndarray): Array of BOLD signal data with shape (time_points, vertices).
    model_prediction (numpy.ndarray): Array of model prediction data with shape (time_points, vertices).
    correction (str, optional): Type of multiple testing correction.
                                Supported methods: 'bonferroni', 'sidak', 'holm-sidak',
                                'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky'.
                                Default is 'fdr_bh'.
    alpha (float or list of floats, optional): The significance level(s) for the tests. Default is 0.01.
    use_fisher (bool, optional): If True, use Fisher's z-score and p-value instead of Pearson's r and p-value. Default is False.

    Returns:
    vertex_results (numpy.ndarray): Array containing the results of linear regression analysis for each vertex.
                                     The shape of the array is (n_output, n_vertex), 
                                     where n_output = slope, intercept, rvalue, pvalue, stderr, trs
                                     + p_values_corrected for each alpha.

    Note:
    The function checks for NaN values in both bold_signal and model_prediction.
    It also identifies and excludes vertices with identical values or containing NaNs.
    """

    # Import 
    import numpy as np
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    def fisher_z(r):
        """Convert Pearson correlation coefficient to Fisher's z-score."""
        if r <= -1 or r >= 1:
            raise ValueError("Pearson correlation coefficient must be in the range (-1, 1).")
        return 0.5 * np.log((1 + r) / (1 - r))
    
    def fisher_z_p_value(z, n):
        """Calculate the two-tailed p-value from Fisher's z-score."""

        se = 1 / np.sqrt(n - 3)  # Standard error
        z_stat = z / se  # Z-statistic
        return 2 * (1 - stats.norm.cdf(np.abs(z_stat)))  # Two-tailed p-value

    if not isinstance(alpha, list):
        alpha = [alpha]
        
    # Check for NaN values in both bold_signal and model_prediction
    nan_mask = np.isnan(model_prediction).any(axis=0) | np.isnan(bold_signal).any(axis=0)

    # Mask for checking identical values along axis 0 in model_prediction
    identical_values_mask = (model_prediction[:-1] == model_prediction[1:]).all(axis=0)

    # Combining nan_mask and identical_values_mask
    invalid_mask = nan_mask | identical_values_mask
    valid_vertices = np.where(~invalid_mask)[0]

    # Define output size
    num_vertices = bold_signal.shape[1]
    num_base_output = 6                         # Number of outputs per vertex (slope, intercept, rvalue, pvalue, stderr, trs)
    num_output = num_base_output + len(alpha)   # Add the desired corrected p-values

    # Array to store results for each vertex
    vertex_results = np.full((num_output, num_vertices), np.nan)  

    # Store p-values before correction
    p_values = np.full(num_vertices, np.nan)

    for i, vert in enumerate(valid_vertices):
        # Perform linear regression
        result = stats.linregress(x=model_prediction[:, vert],
                                  y=bold_signal[:, vert],
                                  alternative='two-sided')
        
        # If using Fisher's z-score
	    if use_fisher:
	        r_value = result.rvalue
	        n = model_prediction.shape[0]  # Number of time points
	        z_score = fisher_z(r_value)
	        p_value = fisher_z_p_value(z_score, n)
	        p_values[vert] = p_value
	        vertex_results[:, vert] = [result.slope, 
	                                   result.intercept, 
	                                   z_score,  # Use Fisher's z-score
	                                   p_value,   # Use Fisher's p-value
	                                   result.stderr,
	                                   n] + [np.nan]*len(alpha)
	    else:
	        # Use standard Pearson values
	        trs = model_prediction.shape[0]
	        vertex_results[:, vert] = [result.slope, 
	                                   result.intercept, 
	                                   result.rvalue, 
	                                   result.pvalue, 
	                                   result.stderr,
	                                   trs] + [np.nan]*len(alpha)
	        p_values[vert] = result.pvalue

    # Apply multiple testing correction
    if correction and valid_vertices.size > 0:
        for n_alphas, alpha_val in enumerate(alpha):
            p_values_corrected = multipletests(p_values[valid_vertices], 
                                               method=correction,
                                               alpha=alpha_val)[1]
            vertex_results[num_base_output + n_alphas, valid_vertices] = p_values_corrected
    
    return vertex_results

def multipletests_surface(pvals, correction='fdr_tsbh', alpha=0.01):
    """
    Perform multiple testing correction for surface data.

    Parameters:
        pvals (numpy.ndarray): Array of p-values.
        correction (str, optional): Method of multiple testing correction. Default is 'fdr_tsbh'.
        alpha (float or list of float, optional): Significance level(s) for the correction.
            Can be a single float or a list of floats for multiple levels. Default is 0.01.

    Returns:
        numpy.ndarray: Array of corrected p-values corresponding to the input p-values.
    """
    # Import
    import numpy as np
    from statsmodels.stats.multitest import multipletests
    
    # If alpha is not already a list, convert it into a list
    if not isinstance(alpha, list):
        alpha_list = [alpha]
    else:
        alpha_list = alpha

    # Check for NaN values in pvals
    nan_mask = np.isnan(pvals)
    valid_vertices = np.where(~nan_mask)[0]

    # Initialize array to store corrected p-values
    corrected_pvals = np.full((len(alpha_list), pvals.shape[0]), np.nan)

    # Perform correction for each specified alpha level
    for n_alpha, alpha_val in enumerate(alpha_list):
        # Perform multiple testing correction and retrieve corrected p-values
        _, p_values_corrected, _, _ = multipletests(pvals[valid_vertices], method=correction, alpha=alpha_val)
        corrected_pvals[n_alpha, valid_vertices] = p_values_corrected

    return corrected_pvals

def median_subject_template(fns): 
    """
    Averages data from different subjects in the same template space.

    Parameters:
        fns (list): List of filenames to be averaged.

    Returns:
        img : Cifti image of the last subject to be used as source_img.
        data_med : The averaged data.
    """
    import numpy as np
    from surface_utils import load_surface
    
    for n_file, fn in enumerate(fns) : 
        print('adding {} to '.format(fn))
        # Load data
        img, data = load_surface(fn=fn)
    
        # Average without nan
        if n_file == 0:
            data_med = np.copy(data)
        else:
            data_med = np.nanmedian(np.array([data_med, data]), axis=0)
            
    return img, data_med

def make_prf_distribution_df(data, rois, max_ecc, grain):
    """
    Load the PRF TSV file and compute the PRF distribution 
    
    Parameters
    ----------
    data: the PRF TSV file
    rois: list of ROIs (Regions Of Interest)
    max_ecc: maximum eccentricity for the Gaussian mesh
    grain: the granularity you want for the Gaussian mesh
    hot_zone_percent: the percentage to define the hot zone (how much of the denser locations you take)
    ci_confidence_level: the confidence level for the confidence interval

    Returns
    -------
    df_distribution: dataframe to use in distribution plot
    """
    import pandas as pd
    import numpy as np
    df_distribution = pd.DataFrame()
    for j, roi in enumerate(rois) :
        # Make df_distribution
        #-------------------
        # Roi data frame
        df_roi = data.loc[data.roi == roi].reset_index()
        if df_roi.empty:
            print(f"[WARNING] No data for ROI: {roi}")
            continue  # skip this ROI
        
        gauss_z_tot = np.zeros((grain,grain)) 
        for vert in range(len(df_roi)):
            # compute the gaussian mesh
            x, y, gauss_z = gaus_2d(gauss_x=df_roi.prf_x[vert],  
                                gauss_y=df_roi.prf_y[vert], 
                                gauss_sd=df_roi.prf_size[vert], 
                                screen_side=max_ecc*2, 
                                grain=grain)
            
            # addition of pRF and ponderation by loo r2
            gauss_z_tot += gauss_z * df_roi.prf_loo_r2[vert]
            
        # Normalisation 
        gauss_z_tot = (gauss_z_tot-gauss_z_tot.min())/(gauss_z_tot.max()-gauss_z_tot.min())
        
        # create the df
        df_distribution_roi = pd.DataFrame()
        df_distribution_roi['roi'] = [roi] * grain
        df_distribution_roi['x'] = x
        df_distribution_roi['y'] = y
        
        gauss_z_tot_df = pd.DataFrame(gauss_z_tot)
        df_distribution_roi = pd.concat([df_distribution_roi, gauss_z_tot_df], axis=1)
        
        if j == 0: df_distribution = df_distribution_roi
        else: df_distribution = pd.concat([df_distribution, df_distribution_roi])
        
    return df_distribution

def make_prf_barycentre_df(df_distribution, rois, max_ecc, grain, hot_zone_percent=0.01):
    """
    Compute the pRF hot zone barycentre
    
    Parameters
    ----------
    df_distribution: df from make_prf_distribution_df
    rois: list of ROIs (Regions Of Interest)
    max_ecc: maximum eccentricity for the Gaussian mesh
    grain: the granularity you want for the Gaussian mesh
    hot_zone_percent: the percentage to define the hot zone (how much of the denser locations you take)
    ci_confidence_level: the confidence level for the confidence interval
        
    Returns
    -------
    df_barycentre: dataframe filtered to use in barycentre plot
    """
    import pandas as pd
    import numpy as np

    for j, roi in enumerate(rois) :
        # Create DataFrame for the region of interest
        df_roi = df_distribution[df_distribution.roi == roi]
        
        # make the two dimensional mesh for z dimension
        exclude_columns = ['roi', 'hemi', 'x', 'y']
        int_columns = df_roi.columns.difference(exclude_columns)
        gauss_z_tot = df_roi[int_columns].values
        
        # Make df_barycentre
        #-------------------
        # Find the 1% higher
        flattened_array = gauss_z_tot.flatten()
        sorted_values = np.sort(flattened_array)[::-1]
        hot_zone_size = int(hot_zone_percent * len(sorted_values))
    
        # make the 2d hot zone idx
        hot_zone_idx = np.unravel_index(np.argsort(flattened_array, axis=None)[-hot_zone_size:], gauss_z_tot.shape)
        
        # Find the barycentre of the top 1% higher
        barycentre_x = np.mean(hot_zone_idx[1])
        barycentre_y = np.mean(hot_zone_idx[0])
        
        # Calculate confidence intervals using bootstrap 
        # num_samples = len(hot_zone_idx[0])
        # lower_ci_x, upper_ci_x = bootstrap_ci_mean(hot_zone_idx[1], n_bootstrap=1000, ci_level=ci_confidence_level)
        # lower_ci_y, upper_ci_y = bootstrap_ci_mean(hot_zone_idx[0], n_bootstrap=1000, ci_level=ci_confidence_level)

        lower_ci_x, upper_ci_x = np.percentile(hot_zone_idx[1], [2.5]), np.percentile(hot_zone_idx[1], [97.5])
        lower_ci_y, upper_ci_y = np.percentile(hot_zone_idx[0], [2.5]), np.percentile(hot_zone_idx[0], [97.5])

        # Convert positions to the correct reference frame
        scale_factor = max_ecc / (grain / 2)
        barycentre_x, barycentre_y = (barycentre_x * scale_factor) - max_ecc, (barycentre_y * scale_factor) - max_ecc
        lower_ci_x, upper_ci_x = (lower_ci_x * scale_factor) - max_ecc, (upper_ci_x * scale_factor) - max_ecc
        lower_ci_y, upper_ci_y = (lower_ci_y * scale_factor) - max_ecc, (upper_ci_y * scale_factor) - max_ecc
        
        # make the df 
        df_barycentre_roi = pd.DataFrame()
        df_barycentre_roi['roi'] = [roi]
        df_barycentre_roi['barycentre_x'] = [barycentre_x]
        df_barycentre_roi['barycentre_y'] = [barycentre_y]
        df_barycentre_roi['lower_ci_x'] = lower_ci_x
        df_barycentre_roi['upper_ci_x'] = upper_ci_x
        df_barycentre_roi['lower_ci_y'] = lower_ci_y
        df_barycentre_roi['upper_ci_y'] = upper_ci_y
        
        if j == 0: df_barycentre = df_barycentre_roi
        else: df_barycentre = pd.concat([df_barycentre, df_barycentre_roi])
        
    return df_barycentre