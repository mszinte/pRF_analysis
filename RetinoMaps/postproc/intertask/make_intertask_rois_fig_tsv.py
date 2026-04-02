"""
-----------------------------------------------------------------------------------------
make_intertask_rois_fig_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs-based CSS tsv
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
CSS analysis tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
python make_intertask_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_intertask_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-hcp1.6mm 327
python make_intertask_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
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
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from maths_utils import make_prf_distribution_df, weighted_nan_median, weighted_nan_percentile, make_prf_barycentre_df

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
glm_settings_path = os.path.join(base_dir, project_dir, "glm-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, glm_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois_methods = analysis_info['rois_methods']
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['prf_amp_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']
subjects_to_group = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
maps_names_pcm = analysis_info['maps_names_pcm']
maps_names_css_stats = analysis_info['maps_names_css_stats']
ecc_size_num_bins = analysis_info['ecc_size_num_bins']
ecc_pcm_num_bins = analysis_info['ecc_pcm_num_bins']
polar_angle_num_bins = analysis_info['polar_angle_num_bins']
ecc_pcm_max = analysis_info['ecc_pcm_max']
ecc_size_max = analysis_info['ecc_size_max']
distribution_max_ecc = analysis_info['distribution_max_ecc']
distribution_mesh_grain = analysis_info['distribution_mesh_grain']
hot_zone_percent = analysis_info['hot_zone_percent']
group_tasks = analysis_info['task_intertask']
categories_to_plot = analysis_info['categories_to_plot']

rsq2use = 'prf_loo_rsq'
avg_method = 'loo-avg'
prf_task_name = analysis_info['prf_task_names'][0]

for tasks in group_tasks : 
    if 'SacVELoc' in tasks: intertask_group = 'SacVE-PurVE-pRF'
    else : intertask_group = 'Sac-Pur-pRF'
    # Format loop
    for format_, extension in zip(formats, extensions):
        print(format_)
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())
        
            # Individual subject analysis
            if 'group' not in subject:
                print('Subject {} is processed'.format(subject))
                
                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    intertask_group, preproc_prep, filtering,
                    normalization, avg_method, rois_method_format)
                
                tsv_dir = "{}/{}/derivatives/pp_data/{}/{}/intertask/tsv".format(main_dir, project_dir, subject, format_)
                tsv_fn = '{}/{}_{}_intertask-deriv.tsv'.format(
                    tsv_dir, subject, fn_spec)
                
                # Exception if no data for one format (e.g template subject)
                if not os.path.isdir(tsv_dir):
                    print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                    continue
                
                data = pd.read_table(tsv_fn, sep="\t")
    
                # replace nan in pcm by non_computed 
                data.loc[data['all'] != 'vision', 
    ['n_neighbor', 'pcm_median', 'vert_geo_dist_median', 'vert_prf_dist_median']] = data.loc[data['all'] != 'vision', 
    ['n_neighbor', 'pcm_median', 'vert_geo_dist_median', 'vert_prf_dist_median']].fillna('non_computed')
    
                # Threshold data (replace by nan)
                if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
                elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
                data.loc[(data.amplitude < amplitude_threshold[0]) |
                         (data.prf_ecc < ecc_threshold[0]) | (data.prf_ecc > ecc_threshold[1]) |
                         (data.prf_size < size_threshold[0]) | (data.prf_size > size_threshold[1]) |
                         (data.prf_n < n_threshold[0]) | (data.prf_n > n_threshold[1]) | 
                         (data[rsq2use] < rsqr_threshold) |
                         (data[stats_col] > stats_threshold)] = np.nan
    
                # loop over categories
                for categorie_to_plot in categories_to_plot : 
                    fn_spec_cate = "task-{}_{}_{}_{}_{}_{}_{}".format(
                        intertask_group, preproc_prep, filtering,
                        normalization, avg_method, rois_method_format, categorie_to_plot)
                    tsv_category_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv/tsv_{}'.format(
                        main_dir, project_dir, subject, format_, categorie_to_plot)
                    os.makedirs(tsv_category_dir, exist_ok=True)
                    data_categorie = data.loc[data['all'] == categorie_to_plot]
                    
                    # Violins
                    # -------
                    df_violins = data_categorie
                    tsv_violins_fn = "{}/{}_{}_intertask_violins.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_violins_fn))
                    df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)
                    
                    # Parameters median
                    # ------------------
                    for num_roi, roi in enumerate(rois):
                        df_roi = data.loc[(data.roi == roi)]
            
                        df_params_median_roi = pd.DataFrame()
                        df_params_median_roi['roi'] = [roi]
                    
                        df_params_median_roi[f'{rsq2use}_weighted_median'] = weighted_nan_median(df_roi[rsq2use], weights=df_roi[rsq2use])
                        df_params_median_roi[f'{rsq2use}_ci_down'] = weighted_nan_percentile(df_roi[rsq2use], df_roi[rsq2use], 25)
                        df_params_median_roi[f'{rsq2use}_ci_up'] = weighted_nan_percentile(df_roi[rsq2use], df_roi[rsq2use], 75)
                        
                        df_params_median_roi['prf_size_weighted_median'] = weighted_nan_median(df_roi.prf_size, weights=df_roi[rsq2use])
                        df_params_median_roi['prf_size_ci_down'] = weighted_nan_percentile(df_roi.prf_size, df_roi[rsq2use], 25)
                        df_params_median_roi['prf_size_ci_up'] = weighted_nan_percentile(df_roi.prf_size, df_roi[rsq2use], 75)
    
                        df_params_median_roi['prf_ecc_weighted_median'] = weighted_nan_median(df_roi.prf_ecc, weights=df_roi[rsq2use])
                        df_params_median_roi['prf_ecc_ci_down'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi[rsq2use], 25)
                        df_params_median_roi['prf_ecc_ci_up'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi[rsq2use], 75)
                        
                        df_params_median_roi['prf_n_weighted_median'] = weighted_nan_median(df_roi.prf_n, weights=df_roi[rsq2use])
                        df_params_median_roi['prf_n_ci_down'] = weighted_nan_percentile(df_roi.prf_n, df_roi[rsq2use], 25)
                        df_params_median_roi['prf_n_ci_up'] = weighted_nan_percentile(df_roi.prf_n, df_roi[rsq2use], 75)
                        
                        pcm_median_col = df_roi.loc[df_roi.pcm_median != 'non_computed', 'pcm_median'].values.astype(float)
                        prf_loo_r2_col = df_roi.loc[df_roi.pcm_median != 'non_computed', '{}'.format(rsq2use)].values.astype(float)
    
                        
                        df_params_median_roi['pcm_median_weighted_median'] = weighted_nan_median(data=pcm_median_col, weights=prf_loo_r2_col)
                        df_params_median_roi['pcm_median_ci_down'] = weighted_nan_percentile(data=pcm_median_col, weights=prf_loo_r2_col, percentile=25)
                        df_params_median_roi['pcm_median_ci_up'] = weighted_nan_percentile(data=pcm_median_col, weights=prf_loo_r2_col, percentile=75)
                
                        if num_roi == 0: df_params_median = df_params_median_roi
                        else: df_params_median = pd.concat([df_params_median, df_params_median_roi])
            
                    tsv_params_median_fn = "{}/{}_{}_intertask_params-median.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_params_median_fn))
                    df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
    
                    # Ecc.size
                    # --------
                    ecc_bins = np.linspace(0, ecc_size_max[1], ecc_size_num_bins+1) 
                    for num_roi, roi in enumerate(rois):
                        df_roi = data.loc[(data.roi == roi)]
                        df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
                        df_ecc_size_bin = pd.DataFrame()
                        df_ecc_size_bin['roi'] = [roi]*ecc_size_num_bins
                        df_ecc_size_bin['num_bins'] = np.arange(ecc_size_num_bins)  
                        df_ecc_size_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x[rsq2use].values)).values
                        df_ecc_size_bin['prf_size_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_size'].values, x[rsq2use].values)).values
                        df_ecc_size_bin[f'{rsq2use}_bins_median'] = np.array(df_bins[rsq2use].median())
                        df_ecc_size_bin['prf_size_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['prf_size'].values, x[rsq2use].values, 75)).values
                        df_ecc_size_bin['prf_size_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['prf_size'].values, x[rsq2use].values, 25)).values
                        if num_roi == 0: df_ecc_size_bins = df_ecc_size_bin
                        else: df_ecc_size_bins = pd.concat([df_ecc_size_bins, df_ecc_size_bin]) 
                    

                    df_ecc_size = df_ecc_size_bins
                    
                    tsv_ecc_size_fn = "{}/{}_{}_intertask_ecc-size.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_ecc_size_fn))
                    df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)
                    
                    # Ecc.pCM
                    # --------
                    data_pcm = data_categorie.loc[data_categorie.pcm_median != 'non_computed']
                    ecc_bins = np.linspace(0, ecc_size_max[1], ecc_pcm_num_bins+1) 
                    for num_roi, roi in enumerate(rois):
                        df_roi = data_pcm.loc[(data.roi == roi)]
                        df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
                        df_ecc_pcm_bin = pd.DataFrame()
                        df_ecc_pcm_bin['roi'] = [roi]*ecc_pcm_num_bins
                        df_ecc_pcm_bin['num_bins'] = np.arange(ecc_pcm_num_bins)
                        df_ecc_pcm_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x[rsq2use].values)).values
                        df_ecc_pcm_bin['prf_pcm_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['pcm_median'].values, x[rsq2use].values)).values
                        df_ecc_pcm_bin['{}_bins_median'.format(rsq2use)] = np.array(df_bins[rsq2use].median())
                        df_ecc_pcm_bin['prf_pcm_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['pcm_median'].values.astype(float), x[rsq2use].values.astype(float), 75)).values
                        df_ecc_pcm_bin['prf_pcm_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['pcm_median'].values.astype(float), x[rsq2use].values.astype(float), 25)).values
                        if num_roi == 0: df_ecc_pcm_bins = df_ecc_pcm_bin
                        else: df_ecc_pcm_bins = pd.concat([df_ecc_pcm_bins, df_ecc_pcm_bin])  
                
                    df_ecc_pcm = df_ecc_pcm_bins
                    
                    df_ecc_pcm = df_ecc_pcm_bins
            
                    tsv_ecc_pcm_fn = "{}/{}_{}_intertask_ecc-pcm.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
                    df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)
                    
                    # Polar angle
                    # -----------
                    theta_slices = np.linspace(0, 360, polar_angle_num_bins, endpoint=False)
                    data['prf_polar_angle'] = np.mod(np.degrees(np.angle(data.polar_real + 1j * data.polar_imag)), 360) 
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for i, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        for j, roi in enumerate(rois): #
                            df = data.loc[(data.roi==roi) & (data.hemi.isin(hemi_values))]
                            if len(df): 
                                df_bins = df.groupby(pd.cut(df['prf_polar_angle'], bins=polar_angle_num_bins))
                                rsq_sum = df_bins[rsq2use].sum()
                            else: 
                                rsq_sum = [np.nan]*polar_angle_num_bins
                
                            df_polar_angle_bin = pd.DataFrame()
                            df_polar_angle_bin['roi'] = [roi]*(polar_angle_num_bins)
                            df_polar_angle_bin['hemi'] = [hemi]*(polar_angle_num_bins)
                            df_polar_angle_bin['num_bins'] = np.arange((polar_angle_num_bins))
                            df_polar_angle_bin['theta_slices'] = np.array(theta_slices)
                            df_polar_angle_bin['rsq_sum'] = np.array(rsq_sum)                        
                            if j == 0 and i == 0: df_polar_angle_bins = df_polar_angle_bin
                            else: df_polar_angle_bins = pd.concat([df_polar_angle_bins, df_polar_angle_bin])  
                                
                    df_polar_angle = df_polar_angle_bins
            
                    tsv_polar_angle_fn = "{}/{}_{}_intertask_polar-angle.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_polar_angle_fn))
                    df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)
                    
                    # Contralaterality
                    # ----------------         
                    for j, roi in enumerate(rois):
                        df_rh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-R')]
                        df_lh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-L')]
                        try: 
                            contralaterality_prct = (sum(df_rh.loc[df_rh.prf_x < 0][rsq2use]) + \
                                                     sum(df_lh.loc[df_lh.prf_x > 0][rsq2use])) / \
                                                    (sum(df_rh[rsq2use]) + sum(df_lh[rsq2use]))
                        except: 
                            contralaterality_prct = np.nan
                        
                        df_contralaterality_roi = pd.DataFrame()
                        df_contralaterality_roi['roi'] = [roi]
                        df_contralaterality_roi['contralaterality_prct'] = np.array(contralaterality_prct)
                
                        if j == 0: df_contralaterality = df_contralaterality_roi
                        else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_roi]) 
            
                    tsv_contralaterality_fn = "{}/{}_{}_intertask_contralaterality.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_contralaterality_fn))
                    df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
                        
                    # Spatial distribution 
                    # --------------------  
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for i, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        data_hemi = data.loc[data.hemi.isin(hemi_values)]
                        df_distribution_hemi = make_prf_distribution_df(
                            data_hemi, rois, distribution_max_ecc, distribution_mesh_grain, rsq2use)
                
                        df_distribution_hemi['hemi'] = [hemi] * len(df_distribution_hemi)
                        if i == 0: df_distribution = df_distribution_hemi
                        else: df_distribution = pd.concat([df_distribution, df_distribution_hemi])
            
                    tsv_distribution_fn = "{}/{}_{}_intertask_distribution.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_distribution_fn))
                    df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)
                    
                    # Spatial distribution hot zone barycentre
                    # ----------------------------------------
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    
                    for i, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
                        df_barycentre_hemi = make_prf_barycentre_df(
                            df_distribution_hemi, rois, distribution_max_ecc, 
                            distribution_mesh_grain, hot_zone_percent=hot_zone_percent)
                        
                        df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
                        if i == 0: df_barycentre = df_barycentre_hemi
                        else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])
                    
                    tsv_barycentre_fn = "{}/{}_{}_intertask_barycentre.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_barycentre_fn))
                    df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)
                               
            # Group Analysis    
            else :
                print('group')
                # loop over categories
                for categorie_to_plot in categories_to_plot :             
                    for i, subject_to_group in enumerate(subjects_to_group):                        
                        tsv_category_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv_{}'.format(
                            main_dir, project_dir, subject_to_group, format_, categorie_to_plot)

                        fn_spec_cate = "task-{}_{}_{}_{}_{}_{}_{}".format(
                            prf_task_name, preproc_prep, filtering, normalization, avg_method, rois_method_format, categorie_to_plot)
                        
                        # Violins
                        # -------
                        tsv_violins_fn = "{}/{}_{}_intertask_violins.tsv".format(
                            tsv_category_dir, subject_to_group, fn_spec_cate)
                        df_violins_indiv = pd.read_table(tsv_violins_fn, sep="\t")
                        if i == 0: df_violins = df_violins_indiv.copy()
                        else: df_violins = pd.concat([df_violins, df_violins_indiv])
                
                        # Parameters median
                        # ------------------
                        # use df_violins
                
                        # Ecc.size
                        # --------
                        tsv_ecc_size_fn = "{}/{}_{}_intertask_ecc-size.tsv".format(
                            tsv_category_dir, subject_to_group, fn_spec_cate)
                        df_ecc_size_indiv = pd.read_table(tsv_ecc_size_fn, sep="\t")
                        if i == 0: df_ecc_size = df_ecc_size_indiv.copy()
                        else: df_ecc_size = pd.concat([df_ecc_size, df_ecc_size_indiv])
                
                        # Ecc.pCM
                        # -------
                        tsv_ecc_pcm_fn = "{}/{}_{}_intertask_ecc-pcm.tsv".format(
                            tsv_category_dir, subject_to_group, fn_spec_cate)
                        df_ecc_pcm_indiv = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
                        if i == 0: df_ecc_pcm = df_ecc_pcm_indiv.copy()
                        else: df_ecc_pcm = pd.concat([df_ecc_pcm, df_ecc_pcm_indiv])
                
                        # Polar angle
                        # -----------
                        tsv_polar_angle_fn = "{}/{}_{}_intertask_polar-angle.tsv".format(
                            tsv_category_dir, subject_to_group, fn_spec_cate)
                        df_polar_angle_indiv = pd.read_table(tsv_polar_angle_fn, sep="\t")
                        if i == 0: df_polar_angle = df_polar_angle_indiv.copy()
                        else: df_polar_angle = pd.concat([df_polar_angle, df_polar_angle_indiv])
                
                        # Contralaterality
                        # ----------------
                        tsv_contralaterality_fn = "{}/{}_{}_intertask_contralaterality.tsv".format(
                            tsv_category_dir, subject_to_group, fn_spec_cate)
                        df_contralaterality_indiv = pd.read_table(tsv_contralaterality_fn, sep="\t")
                        if i == 0: df_contralaterality = df_contralaterality_indiv.copy()
                        else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_indiv])
                        
                        # Spatial distribution 
                        # -------------------
                        tsv_distribution_fn = "{}/{}_{}_intertask_distribution.tsv".format(
                            tsv_category_dir, subject_to_group, fn_spec_cate)
                        df_distribution_indiv = pd.read_table(tsv_distribution_fn, sep="\t")
                        mesh_indiv = df_distribution_indiv.drop(columns=['roi', 'x', 'y', 'hemi']).values
                        others_columns = df_distribution_indiv[['roi', 'x', 'y', 'hemi']]
                        if i == 0: mesh_group = np.expand_dims(mesh_indiv, axis=0)
                        else: mesh_group = np.vstack((mesh_group, np.expand_dims(mesh_indiv, axis=0)))
                      
                    # Median and saving tsv
                    tsv_category_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv/tsv_{}'.format(
                        main_dir, project_dir, subject, format_, categorie_to_plot)
                    os.makedirs(tsv_category_dir, exist_ok=True)
                                    
                    # Violins
                    # -------
                    df_violins = df_violins # no averaging
                    tsv_violins_fn = "{}/{}_{}_intertask_violins.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_violins_fn))
                    df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)
                
                    # Parameters median
                    # ------------------
                    df_params_median = df_violins
                    colnames = ['prf_loo_r2', 'prf_size', 'prf_ecc', 'prf_n', 'pcm_median']
                    if df_violins.empty:
                        df_params_median = pd.DataFrame(columns=['roi', 'subject', 'prf_loo_r2_weighted_median'] + 
                                                        ['{}_weighted_median'.format(col) for col in colnames] + 
                                                        ['{}_ci_down'.format(col) for col in colnames] + 
                                                        ['{}_ci_up'.format(col) for col in colnames])
                    
                    else:
                        df_params_median_indiv = df_params_median.groupby(['roi', 'subject'])[['prf_loo_r2']].apply(
                            lambda x: weighted_nan_median(x['prf_loo_r2'], df_params_median.loc[x.index, 'prf_loo_r2'])).reset_index(name='prf_loo_r2_weighted_median')
                        
                        for colname in colnames[1:]:
                            df_params_median_indiv['{}_weighted_median'.format(colname)] = df_params_median.groupby(['roi', 'subject'])[[colname, 'prf_loo_r2']].apply(
                                lambda x: weighted_nan_median(x[colname], df_params_median.loc[x.index, 'prf_loo_r2'])).reset_index()[0]
                        df_params_med_median = df_params_median_indiv.groupby(['roi'])[[colname + '_weighted_median' for colname in colnames]].median()
                
                        # compute Ci
                        df_params_median_ci = pd.DataFrame()
                        for colname in colnames:
                            df_params_median_ci['{}_ci_down'.format(colname)] = df_params_median_indiv.groupby(['roi']).apply(
                                lambda x: weighted_nan_percentile(x['{}_weighted_median'.format(colname)], x['prf_loo_r2_weighted_median'], 2.5)) 
                            df_params_median_ci['{}_ci_up'.format(colname)] = df_params_median_indiv.groupby(['roi']).apply(
                                lambda x: weighted_nan_percentile(x['{}_weighted_median'.format(colname)], x['prf_loo_r2_weighted_median'], 97.5)) 
                
                    df_params_median = pd.concat([df_params_med_median, df_params_median_ci], axis=1).reset_index()
                    tsv_params_median_fn = "{}/{}_{}_intertask_params-median.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_params_median_fn))
                    df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
                   
                    # Ecc.size
                    # --------
                    df_ecc_size = df_ecc_size.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
                    tsv_ecc_size_fn = "{}/{}_{}_intertask_ecc-size.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_ecc_size_fn))
                    df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)
                
                    # Ecc.pCM
                    # -------
                    df_ecc_pcm = df_ecc_pcm.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
                    tsv_ecc_pcm_fn = "{}/{}_{}_intertask_ecc-pcm.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
                    df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)
                
                    # Polar angle
                    # -----------
                    df_polar_angle = df_polar_angle.groupby(['roi', 'hemi', 'num_bins'], sort=False).median().reset_index()
                    tsv_polar_angle_fn = "{}/{}_{}_intertask_polar-angle.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_polar_angle_fn))
                    df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)
                
                    # Contralaterality
                    # ----------------
                    df_median_contralaterality = df_contralaterality.groupby(['roi'], sort=False).median().reset_index()                
                    df_ci_down_contralaterality = df_contralaterality.groupby(['roi'], sort=False).quantile(0.025).reset_index().rename(columns={'contralaterality_prct': 'ci_down'})
                    df_ci_up_contralaterality = df_contralaterality.groupby(['roi'], sort=False).quantile(0.975).reset_index().rename(columns={'contralaterality_prct': 'ci_up'})
                    df_contralaterality_group = pd.concat([df_median_contralaterality, df_ci_down_contralaterality['ci_down'], df_ci_up_contralaterality['ci_up']], axis=1)
                    tsv_contralaterality_fn = "{}/{}_{}_intertask_contralaterality.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_contralaterality_fn))
                    df_contralaterality_group.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
                    
                    # Spatial distribution 
                    # -------------------
                    tsv_distribution_fn = "{}/{}_{}_intertask_distribution.tsv".format(
                            tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_distribution_fn))
                    median_mesh = np.median(mesh_group, axis=0)
                    df_distribution = pd.DataFrame(median_mesh)
            
                    # Concatenating non-numeric columns back to the dataframe
                    df_distribution = pd.concat([others_columns, df_distribution], axis=1)
                    df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)

                    # Spatial distribution hot zone barycentre 
                    # ----------------------------------------
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for j, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
                        df_barycentre_hemi = make_prf_barycentre_df(
                            df_distribution_hemi, rois, distribution_max_ecc, 
                            distribution_mesh_grain, hot_zone_percent=hot_zone_percent)
                       
                        df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
                        if j == 0: df_barycentre = df_barycentre_hemi
                        else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])
                       
                    tsv_barycentre_fn = "{}/{}_{}_intertask_barycentre.tsv".format(
                        tsv_category_dir, subject, fn_spec_cate)
                    print('Saving tsv: {}'.format(tsv_barycentre_fn))
                    df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)
                    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))    
