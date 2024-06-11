# Figure imports
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import ipdb
import os
deb = ipdb.set_trace

def compute_plot_data(subject, main_dir, project_dir, format_, rois,
                      amplitude_threshold, ecc_threshold, size_threshold, 
                      rsqr_threshold, stats_threshold, pcm_threshold, n_threshold,
                      max_ecc, num_ecc_size_bins, num_ecc_pcm_bins, num_polar_angle_bins,
                      screen_side, gaussian_mesh_grain, hot_zone_percent, subjects_to_group=None):
    """
    Load and compute the data as function of the plot
    
    Parameters
    ----------
    subject: subject string
    main_dir: data main directory
    project_dir: project name
    format_ : format of data
    rois: ROIs
    amplitude_threshold : prf amplitude threshold tupple
    ecc_threshold: prf eccentricity threshold tupple
    size_threshold: prf size threshold tupple
    rsqr_threshold: prf loo r2 threshold tupple
    stats_threshold: prf stats threshold tupple
    n_threshold: prf stats threshold tupple
    pcm_threshold: prf pCM threshold tupple
    max_ecc: maximum eccentricity
    num_ecc_size_bins: number of bins for eccentricty plot for size relationship
    num_ecc_pcm_bins: number of bins for eccentricty plot for pcm relationship
    num_polar_angle_bins: number of bins for eccentricty plot for polar angle
    screen_side: mesh screen side (square) im dva (e.g. 20 dva from -10 to 10 dva)
    gaussian_mesh_grain: The grain you want for the gaussian mesh
    hot_zone_percent: the percentage to define the hot zone (how much of the denser locations you take)
    subjects: list of subject to group
    
    Returns
    -------
    df_roi_area: dataframe to use in surface area plot
    df_violins: dataframe filtered to use in violins plot
    df_ecc_size: dataframe to use in eccentricity size plot
    df_ecc_pcm: dataframe to use in eccentricity size plot
    df_polar_angle: dataframe to use in polar angle plot
    df_contralaterality: dataframe to use in contralaterality plot
    df_params_avg: dataframe to use in parameters average plot
    df_distribution: dataframe to use in prf distribution plot
    """
    from maths_utils import bootstrap_ci_mean, make_prf_distribution_df, make_prf_barycentre_df
    
    if 'group' in subject:
        
        for i, subject_to_group in enumerate(subjects_to_group):
            tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                main_dir, project_dir, subject_to_group, format_)

            # ROI surface areas 
            # -----------------
            tsv_roi_area_fn = "{}/{}_prf_roi_area.tsv".format(tsv_dir, subject_to_group)
            df_roi_area_indiv = pd.read_table(tsv_roi_area_fn, sep="\t")
            if i == 0: df_roi_area = df_roi_area_indiv.copy()
            else: df_roi_area = pd.concat([df_roi_area, df_roi_area_indiv])

            # Violins
            # -------
            tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject_to_group)
            df_violins_indiv = pd.read_table(tsv_violins_fn, sep="\t")
            if i == 0: df_violins = df_violins_indiv.copy()
            else: df_violins = pd.concat([df_violins, df_violins_indiv])

            # Parameters average
            # ------------------
            # use df_violins

            # Ecc.size
            # --------
            tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject_to_group)
            df_ecc_size_indiv = pd.read_table(tsv_ecc_size_fn, sep="\t")
            if i == 0: df_ecc_size = df_ecc_size_indiv.copy()
            else: df_ecc_size = pd.concat([df_ecc_size, df_ecc_size_indiv])

            # Ecc.pCM
            # -------
            tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject_to_group)
            df_ecc_pcm_indiv = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
            if i == 0: df_ecc_pcm = df_ecc_pcm_indiv.copy()
            else: df_ecc_pcm = pd.concat([df_ecc_pcm, df_ecc_pcm_indiv])

            # Polar angle
            # -----------
            tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject_to_group)
            df_polar_angle_indiv = pd.read_table(tsv_polar_angle_fn, sep="\t")
            if i == 0: df_polar_angle = df_polar_angle_indiv.copy()
            else: df_polar_angle = pd.concat([df_polar_angle, df_polar_angle_indiv])

            # Contralaterality
            # ----------------
            tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject_to_group)
            df_contralaterality_indiv = pd.read_table(tsv_contralaterality_fn, sep="\t")
            if i == 0: df_contralaterality = df_contralaterality_indiv.copy()
            else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_indiv])
            
            # Spatial distribution 
            # -------------------
            tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject_to_group)
            df_distribution_indiv = pd.read_table(tsv_distribution_fn, sep="\t")
            if i == 0: df_distribution = df_distribution_indiv.copy()
            else:
                # Identifying numeric columns
                numeric_columns = df_distribution_indiv.select_dtypes(include='number').columns
                non_numeric_columns = df_distribution_indiv.columns.difference(numeric_columns)
                
                # Performing the average only on numeric columns
                df_distribution[numeric_columns] = (df_distribution[numeric_columns] + df_distribution_indiv[numeric_columns]) / 2
                
                # Concatenating non-numeric columns back to the dataframe
                df_distribution[non_numeric_columns] = df_distribution_indiv[non_numeric_columns]

        # Averaging and saving tsv
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
            main_dir, project_dir, subject, format_)
        os.makedirs(tsv_dir, exist_ok=True)
        
        # ROI surface areas 
        # -----------------
        df_roi_area = df_roi_area.groupby(['roi'], sort=False).mean().reset_index()
        tsv_roi_area_fn = "{}/{}_prf_roi_area.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_roi_area_fn))
        df_roi_area.to_csv(tsv_roi_area_fn, sep="\t", na_rep='NaN', index=False)
        
        # Violins
        # -------
        df_violins = df_violins # no averaging
        tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_violins_fn))
        df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)

        # Parameters average
        # ------------------
        df_params_avg = df_violins

        # compute mean and ci
        def weighted_average(df_groupby, column_data, column_weight):
            return (df_groupby[column_data] * df_groupby[column_weight]).sum() / df_groupby[column_weight].sum()

        colnames = ['prf_loo_r2', 'prf_size', 'prf_ecc', 'prf_n', 'pcm']
        df_params_avg_indiv = df_params_avg.groupby(['roi', 'subject'])[['prf_loo_r2']].apply(
                weighted_average, 'prf_loo_r2', 'prf_loo_r2').reset_index(name='prf_loo_r2_weighted_mean')
        for colname in colnames[1:]:            
            df_params_avg_indiv['{}_weighted_mean'.format(colname)] = df_params_avg.groupby(['roi', 'subject'])[[colname, 'prf_loo_r2']].apply(
                    weighted_average, colname, 'prf_loo_r2').reset_index()[0]
        df_params_avg_mean = df_params_avg_indiv.groupby(['roi'])[[colname + '_weighted_mean' for colname in colnames]].mean()
        df_params_avg_ci = pd.DataFrame()
        for colname in colnames:
            df_params_avg_ci['{}_ci_down'.format(colname)] = df_params_avg_indiv.groupby(['roi'])[['{}_weighted_mean'.format(colname)]].apply(lambda x: np.percentile(x, 2.5))
            df_params_avg_ci['{}_ci_up'.format(colname)] = df_params_avg_indiv.groupby(['roi'])[['{}_weighted_mean'.format(colname)]].apply(lambda x: np.percentile(x, 97.5))
        
        df_params_avg = pd.concat([df_params_avg_mean, df_params_avg_ci], axis=1).reset_index()
        tsv_params_avg_fn = "{}/{}_prf_params_avg.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_params_avg_fn))
        df_params_avg.to_csv(tsv_params_avg_fn, sep="\t", na_rep='NaN', index=False)
        
        # Ecc.size
        # --------
        df_ecc_size = df_ecc_size.groupby(['roi', 'num_bins'], sort=False).mean().reset_index()
        tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_size_fn))
        df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

        # Ecc.pCM
        # -------
        df_ecc_pcm = df_ecc_pcm.groupby(['roi', 'num_bins'], sort=False).mean().reset_index()
        tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
        df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)

        # Polar angle
        # -----------
        df_polar_angle = df_polar_angle.groupby(['roi', 'hemi', 'num_bins'], sort=False).mean().reset_index()
        tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_polar_angle_fn))
        df_polar_angle.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)

        # Contralaterality
        # ----------------
        df_contralaterality = df_contralaterality.groupby(['roi'], sort=False).mean().reset_index()
        tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_contralaterality_fn))
        df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
        
        # Spatial distribution 
        # -------------------
        tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_distribution_fn))
        df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)
        
        # Spatial distribution hot zone barycentre 
        # ----------------------------------------
        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        for j, hemi in enumerate(hemis):
            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
            df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
            df_barycentre_hemi = make_prf_barycentre_df(
                df_distribution_hemi, rois, screen_side, gaussian_mesh_grain, hot_zone_percent=hot_zone_percent, ci_confidence_level=0.95)
            
            df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
            if j == 0: df_barycentre = df_barycentre_hemi
            else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])
            
        tsv_barycentre_fn = "{}/{}_prf_barycentre.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_barycentre_fn))
        df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)

    else:
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
            main_dir, project_dir, subject, format_)
        os.makedirs(tsv_dir, exist_ok=True)
        
        tsv_fn = '{}/{}_css-all_derivatives.tsv'.format(tsv_dir, subject)
        data = pd.read_table(tsv_fn, sep="\t")
        
        # keep a raw data df 
        data_raw = data.copy()
    
        # Threshold data (replace by nan)
        if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
        elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
        data.loc[(data.amplitude < amplitude_threshold) |
                 (data.prf_ecc < ecc_threshold[0]) | (data.prf_ecc > ecc_threshold[1]) |
                 (data.prf_size < size_threshold[0]) | (data.prf_size > size_threshold[1]) | 
                 (data.prf_n < n_threshold[0]) | (data.prf_n > n_threshold[1]) | 
                 (data.pcm < pcm_threshold[0]) | (data.pcm > pcm_threshold[1]) |
                 (data.prf_loo_r2 < rsqr_threshold) |
                 (data[stats_col] > stats_threshold)] = np.nan
        data = data.dropna()

        # ROI surface areas 
        # -----------------
        data_raw['vert_area'] = data_raw['vert_area'] / 100 # in cm2
        df_roi_area = data_raw.groupby(['roi'], sort=False)['vert_area'].sum().reset_index()
    
        # Compute the area of FRD 0.05/0.01 vertex in each roi
        df_roi_area['vert_area_corr_pvalue_5pt'] = np.array(data_raw[data_raw['corr_pvalue_5pt'] < 0.05].groupby(
            ['roi'], sort=False)['vert_area'].sum())
        df_roi_area['ratio_corr_pvalue_5pt'] = df_roi_area['vert_area_corr_pvalue_5pt'] / df_roi_area['vert_area'] 
        df_roi_area['vert_area_corr_pvalue_1pt'] = np.array(data_raw[data_raw['corr_pvalue_1pt'] < 0.01].groupby(
            ['roi'], sort=False)['vert_area'].sum())
        df_roi_area['ratio_corr_pvalue_1pt'] = df_roi_area['vert_area_corr_pvalue_1pt'] / df_roi_area['vert_area']         

        # Violins
        # -------
        df_violins = data

        # Parameters average
        # ------------------
        for num_roi, roi in enumerate(rois):
            df_roi = data.loc[(data.roi == roi)]
            df_params_avg_roi = pd.DataFrame()
            df_params_avg_roi['roi'] = [roi]
            df_params_avg_roi['prf_loo_r2_weighted_mean'], _ = np.average(df_roi.prf_loo_r2, weights=df_roi.prf_loo_r2, axis=0, returned=True)
            df_params_avg_roi['prf_loo_r2_ci_down'] = np.percentile(df_roi.prf_loo_r2, [2.5])
            df_params_avg_roi['prf_loo_r2_ci_up'] = np.percentile(df_roi.prf_loo_r2, [97.5])
            df_params_avg_roi['prf_size_weighted_mean'], _ = np.average(df_roi.prf_size, weights=df_roi.prf_loo_r2, axis=0, returned=True)
            df_params_avg_roi['prf_size_ci_down'] = np.percentile(df_roi.prf_size, [2.5])
            df_params_avg_roi['prf_size_ci_up'] = np.percentile(df_roi.prf_size, [97.5])
            df_params_avg_roi['prf_ecc_weighted_mean'], _ = np.average(df_roi.prf_ecc, weights=df_roi.prf_loo_r2, axis=0, returned=True)
            df_params_avg_roi['prf_ecc_ci_down'] = np.percentile(df_roi.prf_ecc, [2.5])
            df_params_avg_roi['prf_ecc_ci_up'] = np.percentile(df_roi.prf_ecc, [97.5])
            df_params_avg_roi['prf_n_weighted_mean'], _ = np.average(df_roi.prf_n, weights=df_roi.prf_loo_r2, axis=0, returned=True)
            df_params_avg_roi['prf_n_ci_down'] = np.percentile(df_roi.prf_n, [2.5])
            df_params_avg_roi['prf_n_ci_up'] = np.percentile(df_roi.prf_n, [97.5])
            df_params_avg_roi['pcm_weighted_mean'], _ = np.average(df_roi.pcm, weights=df_roi.prf_loo_r2, axis=0, returned=True)
            df_params_avg_roi['pcm_ci_down'] = np.percentile(df_roi.pcm, [2.5])
            df_params_avg_roi['pcm_ci_up'] = np.percentile(df_roi.pcm, [97.5])
    
            if num_roi == 0: df_params_avg = df_params_avg_roi
            else: df_params_avg = pd.concat([df_params_avg, df_params_avg_roi])

        # Ecc.size
        # --------
        ecc_bins = np.concatenate(([0],np.linspace(0.4, 1, num_ecc_size_bins)**2 * max_ecc))
        for num_roi, roi in enumerate(rois):
            df_roi = data.loc[(data.roi == roi)]
            df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
            df_ecc_size_bin = pd.DataFrame()
            df_ecc_size_bin['roi'] = [roi]*num_ecc_size_bins
            df_ecc_size_bin['num_bins'] = np.arange(num_ecc_size_bins)
            df_ecc_size_bin['prf_ecc_bins'] = np.array(df_bins['prf_ecc'].mean())
            df_ecc_size_bin['prf_size_bins_mean'] = np.array(df_bins['prf_size'].mean())
            df_ecc_size_bin['prf_loo_r2_bins_mean'] = np.array(df_bins['prf_loo_r2'].mean())
            ci = df_bins['prf_size'].apply(lambda x: bootstrap_ci_mean(x))
            df_ecc_size_bin['prf_size_bins_ci_upper_bound'] = np.array(ci.apply(lambda x: x[1] if not np.isnan(x[1]) else np.nan))
            df_ecc_size_bin['prf_size_bins_ci_lower_bound'] = np.array(ci.apply(lambda x: x[0] if not np.isnan(x[0]) else np.nan))

            if num_roi == 0: df_ecc_size_bins = df_ecc_size_bin
            else: df_ecc_size_bins = pd.concat([df_ecc_size_bins, df_ecc_size_bin])

        df_ecc_size = df_ecc_size_bins

        # Ecc.pCM
        # --------
        data_pcm = data
        ecc_bins = np.concatenate(([0],np.linspace(0.4, 1, num_ecc_pcm_bins)**2 * max_ecc))

        for num_roi, roi in enumerate(rois):
            df_roi = data_pcm.loc[(data.roi == roi)]
            df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
            df_ecc_pcm_bin = pd.DataFrame()
            df_ecc_pcm_bin['roi'] = [roi]*num_ecc_pcm_bins
            df_ecc_pcm_bin['num_bins'] = np.arange(num_ecc_pcm_bins)
            df_ecc_pcm_bin['prf_ecc_bins'] = np.array(df_bins['prf_ecc'].mean())
            df_ecc_pcm_bin['prf_pcm_bins_mean'] = np.array(df_bins['pcm'].mean())
            df_ecc_pcm_bin['prf_loo_r2_bins_mean'] = np.array(df_bins['prf_loo_r2'].mean())
            ci = df_bins['pcm'].apply(lambda x: bootstrap_ci_mean(x))
            df_ecc_pcm_bin['prf_pcm_bins_ci_upper_bound'] = np.array(ci.apply(lambda x: x[1] if not np.isnan(x[1]) else np.nan))
            df_ecc_pcm_bin['prf_pcm_bins_ci_lower_bound'] = np.array(ci.apply(lambda x: x[0] if not np.isnan(x[0]) else np.nan))
            if num_roi == 0: df_ecc_pcm_bins = df_ecc_pcm_bin
            else: df_ecc_pcm_bins = pd.concat([df_ecc_pcm_bins, df_ecc_pcm_bin])

        df_ecc_pcm = df_ecc_pcm_bins

        # Polar angle
        # -----------
        theta_slices = np.linspace(0, 360, num_polar_angle_bins, endpoint=False)
        data['prf_polar_angle'] = np.mod(np.degrees(np.angle(data.polar_real + 1j * data.polar_imag)), 360) 
        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        for i, hemi in enumerate(hemis):
            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
            for j, roi in enumerate(rois): #
                df = data.loc[(data.roi==roi) & (data.hemi.isin(hemi_values))]
                if len(df): 
                    df_bins = df.groupby(pd.cut(df['prf_polar_angle'], bins=num_polar_angle_bins))
                    loo_r2_sum = df_bins['prf_loo_r2'].sum()
                else: loo_r2_sum = [np.nan]*num_polar_angle_bins

                df_polar_angle_bin = pd.DataFrame()
                df_polar_angle_bin['roi'] = [roi]*(num_polar_angle_bins)
                df_polar_angle_bin['hemi'] = [hemi]*(num_polar_angle_bins)
                df_polar_angle_bin['num_bins'] = np.arange((num_polar_angle_bins))
                df_polar_angle_bin['theta_slices'] = np.array(theta_slices)
                df_polar_angle_bin['loo_r2_sum'] = np.array(loo_r2_sum)
                
                if j == 0 and i == 0: df_polar_angle_bins = df_polar_angle_bin
                else: df_polar_angle_bins = pd.concat([df_polar_angle_bins, df_polar_angle_bin])
                    
        df_polar_angle = df_polar_angle_bins
        
        # Contralaterality
        # ----------------         
        for j, roi in enumerate(rois):
            df_rh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-R')]
            df_lh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-L')]
            try: contralaterality_prct = (sum(df_rh.loc[df_rh.prf_x < 0].prf_loo_r2) + \
                                         sum(df_lh.loc[df_lh.prf_x > 0].prf_loo_r2)) / \
                                        (sum(df_rh.prf_loo_r2) + sum(df_lh.prf_loo_r2))
            except: contralaterality_prct = np.nan
            
            df_contralaterality_roi = pd.DataFrame()
            df_contralaterality_roi['roi'] = [roi]
            df_contralaterality_roi['contralaterality_prct'] = np.array(contralaterality_prct)
    
            if j == 0: df_contralaterality = df_contralaterality_roi
            else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_roi])
            
        # Spatial distribution 
        # --------------------  
        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        for i, hemi in enumerate(hemis):
            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
            data_hemi = data.loc[data.hemi.isin(hemi_values)]
            df_distribution_hemi = make_prf_distribution_df(
                data_hemi, rois, screen_side, gaussian_mesh_grain)

            df_distribution_hemi['hemi'] = [hemi] * len(df_distribution_hemi)
            if i == 0: df_distribution = df_distribution_hemi
            else: df_distribution = pd.concat([df_distribution, df_distribution_hemi])
        
        # Spatial distribution hot zone barycentre 
        # ----------------------------------------
        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        for i, hemi in enumerate(hemis):
            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
            df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
            df_barycentre_hemi = make_prf_barycentre_df(
                df_distribution_hemi, rois, screen_side, gaussian_mesh_grain, hot_zone_percent=hot_zone_percent, ci_confidence_level=0.95)
            
            df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
            if i == 0: df_barycentre = df_barycentre_hemi
            else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])
        
        # Saving tsv
        tsv_roi_area_fn = "{}/{}_prf_roi_area.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_roi_area_fn))
        df_roi_area.to_csv(tsv_roi_area_fn, sep="\t", na_rep='NaN', index=False)

        tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_violins_fn))
        df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)

        tsv_params_avg_fn = "{}/{}_prf_params_avg.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_params_avg_fn))
        df_params_avg.to_csv(tsv_params_avg_fn, sep="\t", na_rep='NaN', index=False)

        tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_size_fn))
        df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

        tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
        df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)

        tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_polar_angle_fn))
        df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)

        tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_contralaterality_fn))
        df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
        
        tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_contralaterality_fn))
        df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)
        
        tsv_barycentre_fn = "{}/{}_prf_barycentre.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_barycentre_fn))
        df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)

    return df_roi_area, df_violins, df_ecc_size, df_ecc_pcm, df_polar_angle, df_contralaterality, df_params_avg, df_distribution, df_barycentre
        
def plotly_template(template_specs):
    """
    Define the template for plotly
    Parameters
    ----------
    template_specs : dict
        dictionary contain specific figure settings
    
    Returns
    -------
    fig_template : plotly.graph_objs.layout._template.Template
        Template for plotly figure
    """
    import plotly.graph_objects as go
    fig_template=go.layout.Template()

    # Violin plots
    fig_template.data.violin = [go.Violin(
                                    box_visible=False,
                                    points=False,
                                    # opacity=1,
                                    line_color= "rgba(0, 0, 0, 1)",
                                    line_width=template_specs['plot_width'],
                                    width=0.8,
                                    #marker_symbol='x',
                                    #marker_opacity=1,
                                    hoveron='violins',
                                    meanline_visible=False,
                                    # meanline_color="rgba(0, 0, 0, 1)",
                                    # meanline_width=template_specs['plot_width'],
                                    showlegend=False,
                                    )]

    # Barpolar
    fig_template.data.barpolar = [go.Barpolar(
                                    marker_line_color="rgba(0,0,0,1)",
                                    marker_line_width=template_specs['plot_width'], 
                                    showlegend=False, 
                                    )]
    # Pie plots
    fig_template.data.pie = [go.Pie(textposition=["inside","none"],
                                    # marker_line_color=['rgba(0,0,0,1)','rgba(255,255,255,0)'],
                                    marker_line_width=0,#[template_specs['plot_width'],0],
                                    rotation=0,
                                    direction="clockwise",
                                    hole=0.4,
                                    sort=False,
                                    )]

    # Layout
    fig_template.layout = (go.Layout(# general
                                    font_family=template_specs['font'],
                                    font_size=template_specs['axes_font_size'],
                                    plot_bgcolor=template_specs['bg_col'],

                                    # # x axis
                                    xaxis_visible=True,
                                    xaxis_linewidth=template_specs['axes_width'],
                                    xaxis_color= template_specs['axes_color'],
                                    xaxis_showgrid=False,
                                    xaxis_ticks="outside",
                                    xaxis_ticklen=8,
                                    xaxis_tickwidth = template_specs['axes_width'],
                                    xaxis_title_font_family=template_specs['font'],
                                    xaxis_title_font_size=template_specs['title_font_size'],
                                    xaxis_tickfont_family=template_specs['font'],
                                    xaxis_tickfont_size=template_specs['axes_font_size'],
                                    xaxis_zeroline=False,
                                    xaxis_zerolinecolor=template_specs['axes_color'],
                                    xaxis_zerolinewidth=template_specs['axes_width'],
                                    # xaxis_range=[0,1],
                                    xaxis_hoverformat = '.1f',
                                    
                                    # y axis
                                    yaxis_visible=True,
                                    yaxis_linewidth=template_specs['axes_width'],
                                    yaxis_color= template_specs['axes_color'],
                                    yaxis_showgrid=False,
                                    yaxis_ticks="outside",
                                    yaxis_ticklen=8,
                                    yaxis_tickwidth = template_specs['axes_width'],
                                    yaxis_tickfont_family=template_specs['font'],
                                    yaxis_tickfont_size=template_specs['axes_font_size'],
                                    yaxis_title_font_family=template_specs['font'],
                                    yaxis_title_font_size=template_specs['title_font_size'],
                                    yaxis_zeroline=False,
                                    yaxis_zerolinecolor=template_specs['axes_color'],
                                    yaxis_zerolinewidth=template_specs['axes_width'],
                                    yaxis_hoverformat = '.1f',

                                    # bar polar
                                    polar_radialaxis_visible = False,
                                    polar_radialaxis_showticklabels=False,
                                    polar_radialaxis_ticks='',
                                    polar_angularaxis_visible = False,
                                    polar_angularaxis_showticklabels = False,
                                    polar_angularaxis_ticks = ''
                                    ))

    # Annotations
    fig_template.layout.annotationdefaults = go.layout.Annotation(
                                    font_color=template_specs['axes_color'],
                                    font_family=template_specs['font'],
                                    font_size=template_specs['title_font_size'])

    return fig_template


def prf_roi_area(df_roi_area, fig_width, fig_height, roi_colors):
    """
    Make bar plots of each roi area and the corresponding significative area of pRF  
    
    Parameters
    ----------
    df_roi_area : dataframe for corresponding plot
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    roi_colors : list of rgb colors for plotly
    
    Returns
    -------
    fig : bar plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    
    # General settings
    fig = make_subplots(rows=1, 
                        cols=2, 
                        subplot_titles=['FDR threshold = 0.05', 'FDR threshold = 0.01'],
                       )
    
    # FDR 0.01 
    # All vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area, 
                         text=(df_roi_area['ratio_corr_pvalue_5pt']*100).astype(int).astype(str) + '%',
                         textposition='outside',
                         textangle=-60,
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.2)),
                 row=1, col=1)
 
    # Significant vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area_corr_pvalue_5pt, 
                         showlegend=False, 
                         marker=dict(color=roi_colors)),
                 row=1, col=1)
    
    
    # FDR 0.01 
    # All vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area, 
                         text=(df_roi_area['ratio_corr_pvalue_1pt']*100).astype(int).astype(str) + '%',
                         textposition='outside',
                         textangle=-60,
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.1)),
                 row=1, col=2)
    
    # Significant vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area_corr_pvalue_1pt,
                         showlegend=False, 
                         marker=dict(color=roi_colors)),
                 row=1, col=2)

    # Define parameters
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'))      
    
    fig.update_yaxes(range=[0,100], 
                     showline=True, 
                     nticks=10, 
                     title_text='Surface area (cm<sup>2</sup>)',secondary_y=False)
    
    fig.update_layout(barmode='overlay',
                      height=fig_height, 
                      width=fig_width, 
                      template=fig_template,
                      margin_l=100, 
                      margin_r=50, 
                      margin_t=50, 
                      margin_b=50,
                     )

    # Return outputs
    return fig

def prf_violins_plot(df_violins, fig_width, fig_height, rois, roi_colors):
    """
    Make violins plots for pRF loo_r2, size, n and pcm

    Parameters
    ----------
    df_violins : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    
    Returns
    -------
    fig : violins plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    rows, cols = 2,2
    fig = make_subplots(rows=rows, 
                        cols=cols, 
                        print_grid=False, 
                        vertical_spacing=0.08, 
                        horizontal_spacing=0.05)


    for j, roi in enumerate(rois):
        
        df = df_violins.loc[(df_violins.roi == roi)]
        
        # pRF loo r2
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.prf_loo_r2, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='loo', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 1],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=1, col=1)
                
        # pRF size
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.prf_size, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='size', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 20],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=1, col=2)
        
        # # pRF n
        # fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
        #                         y=df.prf_n, 
        #                         name=roi, 
        #                         opacity=1,
        #                         showlegend=False, 
        #                         legendgroup='n', 
        #                         points=False, 
        #                         scalemode='width', 
        #                         fillcolor=roi_colors[j],
        #                         line_color=roi_colors[j]), 
        #               row=2, col=1)

        # pRF ecc
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.prf_ecc, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='n', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 20],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=2, col=1)
        
        # pcm
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.pcm, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='pcm', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 20],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=2, col=2)
        
        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(showline=True, 
                         range=[0, 1],
                         nticks=10, 
                         title_text='pRF LOO R<sup>2</sup>',
                         row=1, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF size (dva)', 
                         row=1, col=2)
        
        # fig.update_yaxes(showline=True, 
        #                  range=[0, 2], 
        #                  nticks=5, 
        #                  title_text='pRF n', 
        #                  row=2, col=1)

        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF eccentricity (dva)', 
                         row=2, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20],
                         nticks=10, 
                         title_text='pRF pCM (mm/dva)', 
                         row=2, col=2)
        
        fig.update_xaxes(showline=True, 
                         ticklen=0, 
                         linecolor=('rgba(255,255,255,0)'))
        
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=False,
                      legend=dict(orientation="h", 
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1, 
                                  yanchor='top', 
                                  xanchor='left', 
                                  traceorder='normal', 
                                  itemwidth=30), 
                      template=fig_template,
                      margin_l=100, 
                      margin_r=50, 
                      margin_t=100, 
                      margin_b=100)

    return fig

def prf_params_avg_plot(df_params_avg, fig_width, fig_height, rois, roi_colors):
    """
    Make parameters average plots for pRF loo_r2, size, n and pcm

    Parameters
    ----------
    df_params_avg : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    
    Returns
    -------
    fig : parameters average plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    rows, cols = 2,2
    fig = make_subplots(rows=rows, 
                        cols=cols, 
                        print_grid=False, 
                        vertical_spacing=0.08, 
                        horizontal_spacing=0.05)
    
    for j, roi in enumerate(rois):
        
        df = df_params_avg.loc[(df_params_avg.roi == roi)]

        weighted_mean = df.prf_loo_r2_weighted_mean
        ci_up = df.prf_loo_r2_ci_up
        ci_down = df.prf_loo_r2_ci_down
        

        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_mean),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_mean], 
                                              arrayminus=[weighted_mean-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='loo',
                                 showlegend=False), 
                          row=1, col=1)
        
        # pRF size
        weighted_mean = df.prf_size_weighted_mean
        ci_up = df.prf_size_ci_up
        ci_down = df.prf_size_ci_down
        
        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_mean),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_mean], 
                                              arrayminus=[weighted_mean-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='size',
                                 showlegend=False), 
                          row=1, col=2)
                
        # # pRF n
        # weighted_mean = df.prf_n_weighted_mean
        # ci_up = df.prf_n_ci_up
        # ci_down = df.prf_n_ci_down
        
        # fig.add_trace(go.Scatter(x=[roi],
        #                          y=tuple(weighted_mean),
        #                          mode='markers', 
        #                          name=roi,
        #                          error_y=dict(type='data', 
        #                                       array=[ci_up-weighted_mean], 
        #                                       arrayminus=[weighted_mean-ci_down],
        #                                       visible=True, 
        #                                       thickness=3,
        #                                       width=0, 
        #                                       color=roi_colors[j]),
        #                          marker=dict(symbol="square",
        #                                      color=roi_colors[j],
        #                                      size=12, 
        #                                      line=dict(color=roi_colors[j], 
        #                                                width=3)),
        #                          legendgroup='n',
        #                          showlegend=False), 
        #                   row=2, col=1)
        
        # pRF ecc
        weighted_mean = df.prf_ecc_weighted_mean
        ci_up = df.prf_ecc_ci_up
        ci_down = df.prf_ecc_ci_down
        
        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_mean),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_mean], 
                                              arrayminus=[weighted_mean-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='ecc',
                                 showlegend=False), 
                          row=2, col=1)
        
        # pcm
        weighted_mean = df.pcm_weighted_mean
        ci_up = df.pcm_ci_up
        ci_down = df.pcm_ci_down
        
        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_mean),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_mean], 
                                              arrayminus=[weighted_mean-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='pcm',
                                 showlegend=False), 
                          row=2, col=2)

        
        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(showline=True, 
                         range=[0, 1],
                         nticks=10, 
                         title_text='pRF LOO R<sup>2</sup>',
                         row=1, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF size (dva)', 
                         row=1, col=2)
        
        # fig.update_yaxes(showline=True, 
        #                  range=[0, 2], 
        #                  nticks=5, 
        #                  title_text='pRF n', 
        #                  row=2, col=1)

        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF eccentricity (dva)', 
                         row=2, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=10, 
                         title_text='pRF pCM (mm/dva)', 
                         row=2, col=2)
        
        fig.update_xaxes(showline=True, 
                         ticklen=0, 
                         linecolor=('rgba(255,255,255,0)'))
        
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      legend=dict(orientation="h", 
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1, 
                                  yanchor='top', 
                                  xanchor='left', 
                                  traceorder='normal', 
                                  itemwidth=30), 
                      template=fig_template,
                      margin_l=100, 
                      margin_r=50, 
                      margin_t=100, 
                      margin_b=100)

    return fig

def prf_ecc_size_plot(df_ecc_size, fig_width, fig_height, rois, roi_colors, plot_groups, max_ecc):
    """
    Make scatter plot for linear relationship between eccentricity and size

    Parameters
    ----------
    df_ecc_size : A data dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    plot_groups : groups of roi to plot together
    max_ecc : maximum eccentricity 
    
    Returns
    -------
    fig : eccentricy as a function of size plot
    """
    
    from maths_utils import weighted_regression

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(plot_groups)
    fig = make_subplots(rows=rows, cols=cols, print_grid=False)
    
    for l, line_label in enumerate(plot_groups):
        for j, roi in enumerate(line_label):
            
            # Parametring colors
            roi_color = roi_colors[j + l * 3]
            roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"
            
            # Get data
            df = df_ecc_size.loc[(df_ecc_size.roi == roi)]
            ecc_mean = np.array(df.prf_ecc_bins)
            size_mean = np.array(df.prf_size_bins_mean)
            r2_mean = np.array(df.prf_loo_r2_bins_mean)
            size_upper_bound = np.array(df.prf_size_bins_ci_upper_bound)
            size_lower_bound = np.array(df.prf_size_bins_ci_lower_bound)
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_mean, size_mean, r2_mean, model='linear')
            slope_upper, intercept_upper = weighted_regression(ecc_mean[np.where(~np.isnan(size_upper_bound))], 
                                                               size_upper_bound[~np.isnan(size_upper_bound)], 
                                                               r2_mean[np.where(~np.isnan(size_upper_bound))], 
                                                               model='linear')
            slope_lower, intercept_lower = weighted_regression(ecc_mean[np.where(~np.isnan(size_lower_bound))], 
                                                               size_lower_bound[~np.isnan(size_lower_bound)], 
                                                               r2_mean[np.where(~np.isnan(size_lower_bound))], 
                                                               model='linear')

            line_x = np.linspace(ecc_mean[0], ecc_mean[-1], 50)
            line = slope * line_x + intercept
            line_upper = slope_upper * line_x + intercept_upper
            line_lower = slope_lower * line_x + intercept_lower

            fig.add_trace(go.Scatter(x=line_x, y=line, mode='lines', name=roi, legendgroup=roi, 
                                      line=dict(color=roi_color, width=3), showlegend=False), 
                          row=1, col=l+1)

            # Error area
            fig.add_trace(go.Scatter(x=np.concatenate([line_x, line_x[::-1]]), 
                                      y=np.concatenate([list(line_upper), list(line_lower[::-1])]), 
                                      mode='lines', fill='toself', fillcolor=roi_color_opac, 
                                      line=dict(color=roi_color_opac, width=0), showlegend=False), 
                          row=1, col=l+1)

            # Markers
            fig.add_trace(go.Scatter(x=ecc_mean, 
                                     y=size_mean, mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=size_upper_bound - size_mean, 
                                                  arrayminus=size_mean - size_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                      marker=dict(color=roi_color,
                                                  symbol='square',
                                                  size=8, 
                                                  line=dict(color=roi_color, 
                                                            width=3)), 
                                      showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=max_ecc-j*1.5, text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(title_text='pRF size (dva)', row=1, col=1)
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[0, max_ecc], showline=True, row=1, col=l+1)
        fig.update_yaxes(range=[0, max_ecc], showline=True)
        fig.update_layout(height=fig_height, width=fig_width, showlegend=False, template=fig_template,
                         margin_l=100, margin_r=50, margin_t=50, margin_b=100)
        
    return fig

def prf_ecc_pcm_plot(df_ecc_pcm, fig_width, fig_height, rois, roi_colors, plot_groups, max_ecc):
    """
    Make scatter plot for relationship between eccentricity and pCM

    Parameters
    ----------
    df_ecc_pcm : dataframe for the plot
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    plot_groups : groups of roi to plot together
    max_ecc : maximum eccentricity
    
    Returns
    -------
    fig : eccentricy as a function of pcm plot
    """

    from maths_utils import weighted_regression
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(plot_groups)
    fig = make_subplots(rows=rows, cols=cols, print_grid=False)
    
    for l, line_label in enumerate(plot_groups):
        for j, roi in enumerate(line_label):

            # Parametring colors
            roi_color = roi_colors[j + l * 3]
            roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"
            
            # Get data
            df = df_ecc_pcm.loc[(df_ecc_pcm.roi == roi)]
            ecc_mean = np.array(df.prf_ecc_bins)
            pcm_mean = np.array(df.prf_pcm_bins_mean)
            r2_mean = np.array(df.prf_loo_r2_bins_mean)
            pcm_upper_bound = np.array(df.prf_pcm_bins_ci_upper_bound)
            pcm_lower_bound = np.array(df.prf_pcm_bins_ci_lower_bound)
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_mean, pcm_mean, r2_mean, model='pcm')
            
            slope_upper, intercept_upper = weighted_regression(ecc_mean[~np.isnan(pcm_upper_bound)], 
                                                               pcm_upper_bound[~np.isnan(pcm_upper_bound)], 
                                                               r2_mean[~np.isnan(pcm_upper_bound)], 
                                                               model='pcm')
            
            slope_lower, intercept_lower = weighted_regression(ecc_mean[~np.isnan(pcm_lower_bound)], 
                                                               pcm_lower_bound[~np.isnan(pcm_lower_bound)], 
                                                               r2_mean[~np.isnan(pcm_lower_bound)], 
                                                               model='pcm')

            line_x = np.linspace(ecc_mean[0], ecc_mean[-1], 50)
            line = 1 / (slope * line_x + intercept)
            line_upper = 1 / (slope_upper * line_x + intercept_upper)
            line_lower = 1 / (slope_lower * line_x + intercept_lower)

            fig.add_trace(go.Scatter(x=line_x, 
                                     y=line, 
                                     mode='lines', 
                                     name=roi, 
                                     legendgroup=roi, 
                                     line=dict(color=roi_color, width=3), 
                                     showlegend=False), 
                          row=1, col=l+1)

            # Error area
            fig.add_trace(go.Scatter(x=np.concatenate([line_x, line_x[::-1]]),
                                      y=np.concatenate([list(line_upper), list(line_lower[::-1])]), 
                                      mode='lines', fill='toself', fillcolor=roi_color_opac, 
                                      line=dict(color=roi_color_opac, width=0), showlegend=False), 
                          row=1, col=l+1)

            # Markers
            fig.add_trace(go.Scatter(x=ecc_mean, 
                                     y=pcm_mean, 
                                     mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=pcm_upper_bound - pcm_mean, 
                                                  arrayminus=pcm_mean - pcm_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                     marker=dict(color=roi_color, 
                                                 symbol='square',
                                                 size=8, line=dict(color=roi_color,
                                                                   width=3)), 
                                     showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=max_ecc-j*1.5, text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(title_text='pRF cortical magn. (mm/dva)', row=1, col=1)
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[0, max_ecc], showline=True, row=1, col=l+1)
        fig.update_yaxes(range=[0, max_ecc], showline=True)
        fig.update_layout(height=fig_height, width=fig_width, showlegend=False, template=fig_template,
                         margin_l=100, margin_r=50, margin_t=50, margin_b=100)
        
    return fig

def prf_polar_angle_plot(df_polar_angle, fig_width, fig_height, rois, roi_colors, num_polar_angle_bins) :    
    """
    Make polar angle distribution plots
    
    Parameters
    ----------
    df_polar_angle : polar angle dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    num_bins : bins for the polar angle 
     
    Returns
    -------
    figs : a list of three figures
    hemispheres : a list of corresponding hemispheres
    """
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(rois)
    specs = [[{'type': 'polar'}] * cols]
    
    figs = []
    hemispheres = []
    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
    for i, hemi in enumerate(hemis):
        fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs)
            
        for j, roi in enumerate(rois):
            if j == 0: showlegend = True
            else: showlegend = False
    
            # Parts of polar angles and number of voxels in each part
            df = df_polar_angle.loc[(df_polar_angle.roi==roi) & (df_polar_angle.hemi==hemi)]
            
            # barpolar
            fig.add_trace(go.Barpolar(r=df.loo_r2_sum, 
                                      theta=df.theta_slices, 
                                      marker_color=roi_colors[j], 
                                      width=360/(num_polar_angle_bins),
                                      marker_line_color='white', 
                                      marker_line_width=3, 
                                      opacity=1,
                                      showlegend=True,
                                      name=roi, 
                                     ), 
                          row=1, col=j+1)
    
        # Define parameters
        fig.update_polars(angularaxis=dict(visible=False), 
                          radialaxis=dict(visible=False))
        
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          legend=dict(orientation="h", 
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1, 
                                  yanchor='top', 
                                  xanchor='left', 
                                  traceorder='normal', 
                                  itemwidth=30), 
                          template=fig_template,
                          margin_l=50, 
                          margin_r=50, 
                          margin_t=50, 
                          margin_b=50)
                          
        figs.append(fig)
        hemispheres.append(hemi)
        
    return figs, hemispheres

def prf_contralaterality_plot(df_contralaterality, fig_height, fig_width, rois, roi_colors):
    """
    Make contralaterality pie plot
    
    Parameters
    ----------
    df_contralaterality : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
     
    Returns
    -------
    fig : contralaterality figure
    """

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(rois)
    specs = [[{'type': 'pie'}] * cols]    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs)
    
    for j, roi in enumerate(rois):

        df = df_contralaterality.loc[df_contralaterality.roi==roi]
        percentage_total = np.array(df.contralaterality_prct)
        percentage_rest = 1 - percentage_total
        percentage_total = percentage_total.tolist()
        percentage_rest = percentage_rest.tolist()
        values = [percentage_total[0], percentage_rest[0]]


        fig.add_trace(go.Pie(values=values,
                             marker=dict(colors=[roi_colors[j], 'white'])
                            ),
                      row=1, col=j+1)

    # Define parameters
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=False,
                      template=fig_template,
                      margin_l=50, 
                      margin_r=50, 
                      margin_t=50, 
                      margin_b=50)
    
    return fig 

def prf_distribution_plot(df_distribution, fig_height, fig_width, rois, roi_colors, screen_side):
    """
    Make prf distribution contour plot
    
    Parameters
    ----------
    df_distribution : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    screen_side: mesh screen side (square) im dva (e.g. 20 dva from -10 to 10 dva)
     
    Returns
    -------
    fig : distribution figure
    """
    
    # Template settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    fig_template = plotly_template(template_specs)
    
    # General figure settings
    rows, cols = 1, len(rois)
    line_width = 1
    contour_width = 0.5
    
    figs = []
    hemispheres = []
    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
    for i, hemi in enumerate(hemis):  
        fig = make_subplots(rows=rows ,cols=cols)
        for j, roi in enumerate(rois) :
            # Make df roi
            df_roi = df_distribution.loc[(df_distribution.roi == roi) & (df_distribution.hemi == hemi)]

            # make the two dimensional mesh for z dimension
            exclude_columns = ['roi', 'hemi', 'x', 'y']
            int_columns = df_roi.columns.difference(exclude_columns)
            gauss_z_tot = df_roi[int_columns].values
            
            # Contour plot
            fig.add_trace(go.Contour(x=df_roi.x, 
                                     y=df_roi.y, 
                                     z=gauss_z_tot, 
                                     colorscale=[[0, 'white'], [1, roi_colors[j]]],  
                                     showscale=False,  
                                     line=dict(color='black', width=contour_width),  
                                     ncontours=6, 
                                     contours=dict(coloring='fill')
                                    ),row=1, col=j+1)
            
            # x line
            fig.add_trace(go.Scatter(x=[0,0],
                                     y=[-20,20],
                                     mode='lines',
                                     line=dict(dash='2px',color='black', width=line_width)
                                    ),row=1, col=j+1)
            # y line
            fig.add_trace(go.Scatter(x=[-20,20], 
                                     y=[0,0], 
                                     mode='lines', 
                                     line=dict(dash='2px',color='black', width=line_width)),row=1, col=j+1)
            
            # square
            fig.add_shape(type="rect", 
                          x0=-10, 
                          y0=-10, 
                          x1=10, 
                          y1=10, 
                          line=dict(dash='2px',color='black', width=line_width),row=1, col=j+1)
            
        fig.update_xaxes(range=[-20,20], color= ('rgba(255,255,255,0)'))
        fig.update_yaxes(range=[-20,20], color= ('rgba(255,255,255,0)'))
        
        # Define parameters
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          showlegend=False,
                          template=fig_template,
                          margin_l=10, 
                          margin_r=10, 
                          margin_t=100, 
                          margin_b=100)
        figs.append(fig)
        hemispheres.append(hemi)
    return figs, hemispheres

def prf_barycentre_plot(df_barycentre, fig_height, fig_width, rois, roi_colors, screen_side):
    """
    Make prf barycentre plot
    
    Parameters
    ----------
    df_barycentre : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    screen_side: mesh screen side (square) im dva (e.g. 20 dva from -10 to 10 dva)
     
    Returns
    -------
    fig : barycentre figure
    """
    
    # Template settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    fig_template = plotly_template(template_specs)
    
    # General figure settings
    line_width = 1
    fig = go.Figure()
    hemis = ['hemi-L', 'hemi-R']
    for i, hemi in enumerate(hemis): 
        if hemi=='hemi-L': symbol, showlegend = 'square' , True
        elif hemi=='hemi-R': symbol, showlegend = 'circle' , False
        for j, roi in enumerate(rois) :
            # Make df roi
            df_roi = df_barycentre.loc[(df_barycentre.roi == roi) & (df_barycentre.hemi == hemi)]    
    
            # barycentre position
            fig.add_trace(go.Scatter(x=df_roi.barycentre_x, 
                                     y=df_roi.barycentre_y, 
                                     mode='markers', 
                                     name = roi,
                                     marker=dict(symbol=symbol, 
                                                 color=roi_colors[j], 
                                                 size=12),
                                     error_x=dict(type='data', 
                                                  array=[df_roi.upper_ci_x - df_roi.barycentre_x], 
                                                  arrayminus=[df_roi.barycentre_x - df_roi.lower_ci_x],
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_colors[j]),
                                     error_y=dict(type='data', 
                                                  array=[df_roi.upper_ci_y - df_roi.barycentre_y], 
                                                  arrayminus=[df_roi.barycentre_y - df_roi.lower_ci_y],
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_colors[j]),
                                 showlegend=showlegend))
        # Center lignes
        fig.add_trace(go.Scatter(x=[0,0], 
                                 y=[-screen_side, screen_side], 
                                 mode='lines', 
                                 showlegend=False, 
                                 line=dict(dash='2px',color='grey', width=line_width)))
        
        fig.add_trace(go.Scatter(x=[-screen_side,screen_side], 
                                 y=[0,0], 
                                 mode='lines', 
                                 showlegend=False,
                                 line=dict(dash='2px',color='grey', width=line_width)))
        
        # Add squares 
        for position in [2,4,6,8,10]:
            fig.add_shape(type="rect", 
                          x0=-position, 
                          y0=-position, 
                          x1=position, 
                          y1=position, 
                          line=dict(dash='2px',color='grey', width=line_width))
        # Add annotations 
        fig.add_trace(go.Scatter(x=[0, 0, 0, 0, 0], 
                                 y=[2.2, 4.2, 6.2, 8.2, 10.2], 
                                 showlegend=False, 
                                 text=["2 dva", 
                                       "4 dva", 
                                       "6 dva",
                                       "8 dva", 
                                       "10 dva"], 
                                 mode="text", 
                                 textfont=dict(size=10)))
    
    fig.update_yaxes(range=[-12,12],color= ('rgba(255,255,255,0)'))
    fig.update_xaxes(range=[-12,12],color= ('rgba(255,255,255,0)'))

    # Define parameters
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=True,
                      template=fig_template,
                      margin_l=570, 
                      margin_r=570, 
                      margin_t=50, 
                      margin_b=50)
        
    return fig

def categories_proportions_roi_plot(data, subject, fig_height, fig_width):
    data = data.copy()
    filtered_data = data[data['stats_final'] != 'non_responding']
    
    # Sort categories
    categories_order = ['vision', 'vision_and_pursuit_and_saccade', 'pursuit_and_saccade', 'vision_and_saccade', 'vision_and_pursuit', 'saccade', 'pursuit']
    filtered_data['stats_final'] = pd.Categorical(filtered_data['stats_final'], categories=categories_order, ordered=True)
    filtered_data = filtered_data.sort_values(['rois', 'stats_final'])

    
    #  Defines colors settings 
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    # To write make the percent visible only for choose categories 
    percent_color =  {'pursuit': 'rgba(255,255,255,0)', 
                      'saccade': 'rgba(255,255,255,0)', 
                      'pursuit_and_saccade': 'rgba(0, 0, 0, 1)', 
                      'vision': 'rgba(0, 0, 0, 1)', 
                      'vision_and_pursuit': 'rgba(255,255,255,0)', 
                      'vision_and_saccade': 'rgba(255,255,255,0)', 
                      'vision_and_pursuit_and_saccade': 'rgba(0, 0, 0, 1)'}
    
    categorie_color_map = {'pursuit': 'rgba(255,255,255,0)', 
                           'saccade': 'rgba(255,255,255,0)', 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': 'rgba(255,255,255,0)', 
                           'vision_and_saccade': 'rgba(255,255,255,0)', 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    rois = pd.unique(data.rois)
    #  Make the subplot
    # fig_height, fig_width = 300, 1920
    rows = 2 
    cols =len(rois)
    specs = [[{'type': 'domain'}] * cols,  [{'type': 'xy'}] * cols]
    
    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs, row_heights=[1,0.05])
    
    
    for i, roi in enumerate(rois):
        df_rois = filtered_data.loc[filtered_data.rois == roi]
        #  Colors for categories 
        categorie_colors = [categorie_color_map[label] for label in df_rois.stats_final]
        #  Colors for the percentages 
        percentage_colors = [percent_color[label] for label in df_rois.stats_final]
        
        
        fig.add_trace(go.Pie(labels=df_rois.stats_final, 
                             values=df_rois.vertex_surf, 
                             showlegend=False, 
                             sort=False,
                             textinfo='percent',
                             textposition='inside',
                             direction='clockwise',
                             name= roi,    
                             marker=dict(colors=categorie_colors),
                             insidetextfont=dict(color=percentage_colors),
                             hole=0.3), 
                      row=1, col=i+1)
        
        fig.add_annotation(text=roi, 
                           yshift =10,
                           showarrow=False, 
                           font=dict(size=13,color=roi_colors[i]), 
                           row=2, col=i+1)
    
    
        
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      template='simple_white')  
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(color='rgba(255,255,255,0)'))
    
    fig.update_yaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(color='rgba(255,255,255,0)'))
    
    return fig

def surface_rois_categories_plot(data, subject, fig_height, fig_width):   
    data = data.copy()
    #  Defines colors settings 
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    
    
    
    categorie_color_map = {'non_responding': stats_categories_colors[0], 
                           'pursuit': stats_categories_colors[1], 
                           'saccade': stats_categories_colors[2], 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': stats_categories_colors[5], 
                           'vision_and_saccade': stats_categories_colors[6], 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    
    #  grpup df 
    group_df_rois = data.groupby(['rois'], sort=False)['vertex_surf'].sum().reset_index()
    group_rois_categories = data.groupby(['rois', 'stats_final'], sort=False)['vertex_surf'].sum().reset_index()
    
    #  Make subplot 
    # fig_height, fig_width = 1080, 1920
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    
    
    fig.add_trace(go.Bar(x=group_df_rois.rois, 
                         y=group_df_rois.vertex_surf, 
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.1)), 
                  secondary_y=False)
    
    #Choose categories to plot
    stats_categories = ['vision',  'vision_and_pursuit_and_saccade','pursuit_and_saccade']
    for categorie in stats_categories:
        df = group_rois_categories.loc[group_rois_categories.stats_final == categorie]
    
        fig.add_trace(go.Bar(x=df.rois, 
                             y=df.vertex_surf, 
                             name=categorie,  
                             legendgroup=categorie, 
                             marker_color=categorie_color_map[categorie]), 
                      secondary_y=True) 
    
    fig.update_layout(yaxis2=dict(overlaying='y',
                                  side='right',
                                  range=[0, 5000],  
                                  showticklabels=False, 
                                  ticklen=0, 
                                  linecolor=('rgba(255,255,255,0)')))
    
    
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(size=12))      
    
    fig.update_yaxes(range=[0,5000], 
                     nticks=5, 
                     title_text='Surface in mm<sup>2</sup>',secondary_y=False)
    
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      barmode='stack',
                      showlegend=True, 
                      template='simple_white')  
    
    
    return fig 


def surface_rois_all_categories_plot(data, subject, fig_height, fig_width):  
    data = data.copy()
    #  Defines colors settings 
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    
    
    
    categorie_color_map = {'non_responding': stats_categories_colors[0], 
                           'pursuit': stats_categories_colors[1], 
                           'saccade': stats_categories_colors[2], 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': stats_categories_colors[5], 
                           'vision_and_saccade': stats_categories_colors[6], 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    #  grpup df 
    group_df = data.groupby(['rois', 'stats_final'], sort=False)['vertex_surf'].sum().reset_index()
    
    #  Figure settings
    # fig_height, fig_width = 1080, 1920
    fig = go.Figure()
    
    
    stats_categories= ['non_responding', 'vision','vision_and_pursuit_and_saccade', 'pursuit_and_saccade', 'pursuit', 'saccade', 'vision_and_pursuit', 'vision_and_saccade']
    for categorie in stats_categories:
        df = group_df.loc[group_df.stats_final == categorie]
    
        fig.add_trace(go.Bar(x=df.rois, 
                             y=df.vertex_surf, 
                             name=categorie,  
                             legendgroup=categorie, 
                             marker_color=categorie_color_map[categorie])) 
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(size=12))      
    
    fig.update_yaxes(range=[0,5000], 
                     nticks=5, 
                     title_text='Surface in mm<sup>2</sup>')
    
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      barmode='stack',
                      showlegend=True, 
                      template='simple_white')  
    
    return fig 