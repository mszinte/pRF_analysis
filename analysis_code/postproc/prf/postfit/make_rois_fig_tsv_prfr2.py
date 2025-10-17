"""
-----------------------------------------------------------------------------------------
make_rois_fig_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs figure specific TSV
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
sys.argv[5]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
CSS analysis tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] 
                        [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/

python make_rois_fig_tsv.py /scratch/mszinte/data MotConf sub-01 327
python make_rois_fig_tsv.py /scratch/mszinte/data MotConf sub-170k 327
python make_rois_fig_tsv.py /scratch/mszinte/data MotConf group 327

python make_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps group 327

python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf sub-01 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf sub-170k 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_patient 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_aniso 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_mixed 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_strab 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_control 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_excluded 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
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
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from maths_utils import make_prf_distribution_df, weighted_nan_median, weighted_nan_percentile

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
if len(sys.argv) > 5: output_folder = sys.argv[5]
else: output_folder = "prf"

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['amplitude_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']
subjects_to_group = analysis_info['subjects']
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']

fig_settings_path = os.path.join(base_dir, project_dir, "figure_settings.json")
with open(fig_settings_path) as f:
    json_s = f.read()
    figure_info = json.loads(json_s)
num_ecc_size_bins = figure_info['num_ecc_size_bins']
num_ecc_pcm_bins = figure_info['num_ecc_pcm_bins']
num_polar_angle_bins = figure_info['num_polar_angle_bins']
max_ecc = figure_info['max_ecc']
screen_side = figure_info['screen_side']
gaussian_mesh_grain = figure_info['gaussian_mesh_grain']
hot_zone_percent = figure_info['hot_zone_percent']

# Format loop
for format_, extension in zip(formats, extensions):
    print(format_)
    
    # Individual subject analysis
    if 'group' not in subject:
        print('Subject {} is processed'.format(subject))
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/{}/tsv'.format(
            main_dir, project_dir, subject, format_, output_folder)
        os.makedirs(tsv_dir, exist_ok=True)
        
        tsv_fn = '{}/{}_css-all_derivatives.tsv'.format(tsv_dir, subject)
        data = pd.read_table(tsv_fn, sep="\t")

       
        # Keep a raw data df 
        data_raw = data.copy()
        
        # Threshold data (replace by nan)
        # if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
        # elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
        # data.loc[(data.amplitude < amplitude_threshold) |
        #          (data.prf_ecc < ecc_threshold[0]) | (data.prf_ecc > ecc_threshold[1]) |
        #          (data.prf_size < size_threshold[0]) | (data.prf_size > size_threshold[1]) | 
        #          (data.prf_n < n_threshold[0]) | (data.prf_n > n_threshold[1]) | 
        #          (data.prf_rsq < rsqr_threshold) |
        #          (data[stats_col] > stats_threshold)] = np.nan
        #data = data.dropna()
        
        
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
        
        tsv_roi_area_fn = "{}/{}_prf_roi_area.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_roi_area_fn))
        df_roi_area.to_csv(tsv_roi_area_fn, sep="\t", na_rep='NaN', index=False)
    
        # Violins
        # -------
        df_violins = data
        
        tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_violins_fn))
        df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)
        
        # Parameters median
        # ------------------
        # for num_roi, roi in enumerate(rois):
        #     df_roi = data.loc[(data.roi == roi)]

        #     df_params_median_roi = pd.DataFrame()
        #     df_params_median_roi['roi'] = [roi]
        
        #     #df_params_median_roi['prf_rsq_weighted_median'] = weighted_nan_median(df_roi.prf_rsq, weights=df_roi.prf_rsq)
        #     df_params_median_roi['prf_rsq_weighted_mediann'] = weighted_nan_median(df_roi.prf_rsq, weights=df_roi.prf_rsq)
        #     df_params_median_roi['prf_rsq_ci_down'] = weighted_nan_percentile(df_roi.prf_rsq, df_roi.prf_rsq, 25)
        #     df_params_median_roi['prf_rsq_ci_up'] = weighted_nan_percentile(df_roi.prf_rsq, df_roi.prf_rsq, 75)
            
        #     df_params_median_roi['prf_size_weighted_median'] = weighted_nan_median(df_roi.prf_size, weights=df_roi.prf_rsq)
        #     df_params_median_roi['prf_size_ci_down'] = weighted_nan_percentile(df_roi.prf_size, df_roi.prf_rsq, 25)
        #     df_params_median_roi['prf_size_ci_up'] = weighted_nan_percentile(df_roi.prf_size, df_roi.prf_rsq, 75)
            
        #     df_params_median_roi['prf_ecc_weighted_median'] = weighted_nan_median(df_roi.prf_ecc, weights=df_roi.prf_rsq)
        #     df_params_median_roi['prf_ecc_ci_down'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi.prf_rsq, 25)
        #     df_params_median_roi['prf_ecc_ci_up'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi.prf_rsq, 75)
            
        #     df_params_median_roi['prf_n_weighted_median'] = weighted_nan_median(df_roi.prf_n, weights=df_roi.prf_rsq)
        #     df_params_median_roi['prf_n_ci_down'] = weighted_nan_percentile(df_roi.prf_n, df_roi.prf_rsq, 25)
        #     df_params_median_roi['prf_n_ci_up'] = weighted_nan_percentile(df_roi.prf_n, df_roi.prf_rsq, 75)
             
        #     df_params_median_roi['pcm_median_weighted_median'] = weighted_nan_median(df_roi.pcm_median, weights=df_roi.prf_rsq)
        #     df_params_median_roi['pcm_median_ci_down'] = weighted_nan_percentile(df_roi.pcm_median, df_roi.prf_rsq, 25)
        #     df_params_median_roi['pcm_median_ci_up'] = weighted_nan_percentile(df_roi.pcm_median, df_roi.prf_rsq, 75)
    
        #     if num_roi == 0: df_params_median = df_params_median_roi
        #     else: df_params_median = pd.concat([df_params_median, df_params_median_roi])

        # tsv_params_median_fn = "{}/{}_prf_params_median.tsv".format(tsv_dir, subject)
        # print('Saving tsv: {}'.format(tsv_params_median_fn))
        # df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
        
        # weight by rsq instead prf_rsq
        for num_roi, roi in enumerate(rois):
            df_roi = data.loc[(data.roi == roi)]

            df_params_median_roi = pd.DataFrame()
            df_params_median_roi['roi'] = [roi]

            # Use prf_rsq for both median and weights
            df_params_median_roi['prf_rsq_weighted_median'] = weighted_nan_median(df_roi.prf_rsq, weights=df_roi.prf_rsq)
            df_params_median_roi['prf_rsq_ci_down'] = weighted_nan_percentile(df_roi.prf_rsq, df_roi.prf_rsq, 25)
            df_params_median_roi['prf_rsq_ci_up'] = weighted_nan_percentile(df_roi.prf_rsq, df_roi.prf_rsq, 75)

            df_params_median_roi['prf_size_weighted_median'] = weighted_nan_median(df_roi.prf_size, weights=df_roi.prf_rsq)
            df_params_median_roi['prf_size_ci_down'] = weighted_nan_percentile(df_roi.prf_size, df_roi.prf_rsq, 25)
            df_params_median_roi['prf_size_ci_up'] = weighted_nan_percentile(df_roi.prf_size, df_roi.prf_rsq, 75)

            df_params_median_roi['prf_ecc_weighted_median'] = weighted_nan_median(df_roi.prf_ecc, weights=df_roi.prf_rsq)
            df_params_median_roi['prf_ecc_ci_down'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi.prf_rsq, 25)
            df_params_median_roi['prf_ecc_ci_up'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi.prf_rsq, 75)

            df_params_median_roi['prf_n_weighted_median'] = weighted_nan_median(df_roi.prf_n, weights=df_roi.prf_rsq)
            df_params_median_roi['prf_n_ci_down'] = weighted_nan_percentile(df_roi.prf_n, df_roi.prf_rsq, 25)
            df_params_median_roi['prf_n_ci_up'] = weighted_nan_percentile(df_roi.prf_n, df_roi.prf_rsq, 75)

            #df_params_median_roi['pcm_median_weighted_median'] = weighted_nan_median(df_roi.pcm_median, weights=df_roi.prf_rsq)
            #df_params_median_roi['pcm_median_ci_down'] = weighted_nan_percentile(df_roi.pcm_median, df_roi.prf_rsq, 25)
            #df_params_median_roi['pcm_median_ci_up'] = weighted_nan_percentile(df_roi.pcm_median, df_roi.prf_rsq, 75)

            if num_roi == 0:
                df_params_median = df_params_median_roi
            else:
                df_params_median = pd.concat([df_params_median, df_params_median_roi])

        tsv_params_median_fn = f"{tsv_dir}/{subject}_prf_params_median.tsv"
        print(f"Saving tsv: {tsv_params_median_fn}")
        df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)

        # Ecc.size
        # --------
        # ecc_bins = np.concatenate(([0],np.linspace(0.26, 1, num_ecc_size_bins)**2 * max_ecc))
        # for num_roi, roi in enumerate(rois):
        #     df_roi = data.loc[(data.roi == roi)]
        #     df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
        #     df_ecc_size_bin = pd.DataFrame()
        #     df_ecc_size_bin['roi'] = [roi]*num_ecc_size_bins
        #     df_ecc_size_bin['num_bins'] = np.arange(num_ecc_size_bins)
        #     df_ecc_size_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x['prf_rsq'].values)).values

        #     df_ecc_size_bin['prf_size_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_size'].values, x['prf_rsq'].values)).values
        #     df_ecc_size_bin['prf_rsq_bins_median'] = np.array(df_bins['prf_rsq'].median())
        #     df_ecc_size_bin['prf_size_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['prf_size'].values, x['prf_rsq'].values, 75)).values
        #     df_ecc_size_bin['prf_size_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['prf_size'].values, x['prf_rsq'].values, 25)).values
        #     if num_roi == 0: df_ecc_size_bins = df_ecc_size_bin
        #     else: df_ecc_size_bins = pd.concat([df_ecc_size_bins, df_ecc_size_bin]) 
        
        # df_ecc_size = df_ecc_size_bins
        
        # tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject)
        # print('Saving tsv: {}'.format(tsv_ecc_size_fn))
        # df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

        ecc_bins = np.concatenate(([0], np.linspace(0.26, 1, num_ecc_size_bins)**2 * max_ecc))

        # Define the R² column once to keep it clean
        r2_col = 'prf_rsq' if 'prf_rsq' in data.columns else 'prf_rsq'

        for num_roi, roi in enumerate(rois):
            df_roi = data.loc[(data.roi == roi)]
            df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
            df_ecc_size_bin = pd.DataFrame()
            df_ecc_size_bin['roi'] = [roi] * num_ecc_size_bins
            df_ecc_size_bin['num_bins'] = np.arange(num_ecc_size_bins)

            # --- Eccentricity per bin (weighted by R²) ---
            df_ecc_size_bin['prf_ecc_bins'] = df_bins.apply(
                lambda x: weighted_nan_median(x['prf_ecc'].values, x[r2_col].values)
            ).values

            # --- Size per bin (weighted by R²) ---
            df_ecc_size_bin['prf_size_bins_median'] = df_bins.apply(
                lambda x: weighted_nan_median(x['prf_size'].values, x[r2_col].values)
            ).values

            # --- R² per bin (median, not weighted) ---
            df_ecc_size_bin['prf_rsq_bins_median'] = np.array(df_bins[r2_col].median())

            # --- Size CI bounds per bin ---
            df_ecc_size_bin['prf_size_bins_ci_upper_bound'] = df_bins.apply(
                lambda x: weighted_nan_percentile(x['prf_size'].values, x[r2_col].values, 75)
            ).values

            df_ecc_size_bin['prf_size_bins_ci_lower_bound'] = df_bins.apply(
                lambda x: weighted_nan_percentile(x['prf_size'].values, x[r2_col].values, 25)
            ).values

            if num_roi == 0:
                df_ecc_size_bins = df_ecc_size_bin
            else:
                df_ecc_size_bins = pd.concat([df_ecc_size_bins, df_ecc_size_bin])

        df_ecc_size = df_ecc_size_bins

        tsv_ecc_size_fn = f"{tsv_dir}/{subject}_prf_ecc_size.tsv"
        print(f"Saving tsv: {tsv_ecc_size_fn}")
        df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

        
        # # Ecc.pCM
        # --------
        # data_pcm = data
        # ecc_bins = np.concatenate(([0],np.linspace(0.26, 1, num_ecc_pcm_bins)**2 * max_ecc))
        # for num_roi, roi in enumerate(rois):
        #     df_roi = data_pcm.loc[(data.roi == roi)]
        #     df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
        #     df_ecc_pcm_bin = pd.DataFrame()
        #     df_ecc_pcm_bin['roi'] = [roi]*num_ecc_pcm_bins
        #     df_ecc_pcm_bin['num_bins'] = np.arange(num_ecc_pcm_bins)
        #     df_ecc_pcm_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x['prf_rsq'].values)).values
        #     df_ecc_pcm_bin['prf_pcm_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['pcm_median'].values, x['prf_rsq'].values)).values
        #     df_ecc_pcm_bin['prf_rsq_bins_median'] = np.array(df_bins['prf_rsq'].median())
        #     df_ecc_pcm_bin['prf_pcm_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['pcm_median'].values, x['prf_rsq'].values, 75)).values
        #     df_ecc_pcm_bin['prf_pcm_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['pcm_median'].values, x['prf_rsq'].values, 25)).values
        #     if num_roi == 0: df_ecc_pcm_bins = df_ecc_pcm_bin
        #     else: df_ecc_pcm_bins = pd.concat([df_ecc_pcm_bins, df_ecc_pcm_bin])  
    
        # df_ecc_pcm = df_ecc_pcm_bins

        # tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject)
        # print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
        # df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)
        
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
                    prf_rsq_sum = df_bins['prf_rsq'].sum()
                else: prf_rsq_sum = [np.nan]*num_polar_angle_bins
    
                df_polar_angle_bin = pd.DataFrame()
                df_polar_angle_bin['roi'] = [roi]*(num_polar_angle_bins)
                df_polar_angle_bin['hemi'] = [hemi]*(num_polar_angle_bins)
                df_polar_angle_bin['num_bins'] = np.arange((num_polar_angle_bins))
                df_polar_angle_bin['theta_slices'] = np.array(theta_slices)
                df_polar_angle_bin['prf_rsq_sum'] = np.array(prf_rsq_sum)
                
                if j == 0 and i == 0: df_polar_angle_bins = df_polar_angle_bin
                else: df_polar_angle_bins = pd.concat([df_polar_angle_bins, df_polar_angle_bin])  
                    
        df_polar_angle = df_polar_angle_bins

        tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_polar_angle_fn))
        df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)

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

        tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_distribution_fn))
        df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)
        
        # Contralaterality
        # ----------------         
        for j, roi in enumerate(rois):
            df_rh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-R')]
            df_lh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-L')]
            try: contralaterality_prct = (sum(df_rh.loc[df_rh.prf_x < 0].prf_rsq) + \
                                         sum(df_lh.loc[df_lh.prf_x > 0].prf_rsq)) / \
                                        (sum(df_rh.prf_rsq) + sum(df_lh.prf_rsq))
            except: contralaterality_prct = np.nan
            
            df_contralaterality_roi = pd.DataFrame()
            df_contralaterality_roi['roi'] = [roi]
            df_contralaterality_roi['contralaterality_prct'] = np.array(contralaterality_prct)
    
            if j == 0: df_contralaterality = df_contralaterality_roi
            else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_roi]) 

        tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_contralaterality_fn))
        df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
            
       
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))          
        