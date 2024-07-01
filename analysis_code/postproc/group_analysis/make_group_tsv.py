"""
-----------------------------------------------------------------------------------------
make_group_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
make a sub-all tsv
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# sub-all tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/group_analysis
2. run python command
python make_sub_all_tsv.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python make_group_tsv.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug 
import ipdb
deb = ipdb.set_trace

# general imports
import os
import sys
import json
import numpy as np
import pandas as pd 

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois_list = analysis_info['rois']


#  Defind threshold 
ecc_th = [0, 30]
size_th= [0,30 ]
rsq_th = [0, 1]
pcm_th = [0, 30]
n_th = [0, 10]




# Iterate over formats and data types
for format_, extension in zip(formats, extensions):
    print(format_)
    if format_ == 'fsnative':
        data_types = ['derivatives']
    elif format_ == '170k':
        data_types = ['derivatives', 'derivatives_group']
        
    # Initialize DataFrames to store aggregated data
    data_all = pd.DataFrame()
    data_mean_all = pd.DataFrame() 
    for data_type in data_types:
        print(data_type)

        
        # Iterate over subjects
        for subject in subjects:
            print('adding {}'.format(subject))
            tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(main_dir, project_dir, subject, format_)

            # Read data from the subject's file
            data = pd.read_table('{}/{}_css-prf_{}.tsv'.format(tsv_dir,subject,data_type))
            
            # Apply thresholds to filter out invalid data
            data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
                     (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) | 
                     (data.prf_n < n_th[0]) | (data.prf_n > n_th[1]) | 
                     (data.prf_loo_r2 <= rsq_th[0])] = np.nan
            data = data.dropna()

            # Concatenate data from all subjects
            data_all = pd.concat([data_all, data], ignore_index=True)

            # Calculate mean data ponderate by r2 for each subject by bins of ecc for both hemi and all categories
            numeric_columns = data.select_dtypes(include='number').columns
            grouped_data = data.groupby(['hemi', 'rois', 'stats_final', pd.cut(data['prf_ecc'], bins=np.arange(0, 20, 1))])
            data_mean_subject = grouped_data.apply(lambda x: (x[numeric_columns].drop(columns='prf_rsq').mul(x['prf_rsq'], axis=0).sum() / x['prf_rsq'].sum())).assign(prf_rsq=grouped_data['prf_rsq'].mean()).reset_index(level=['hemi', 'rois','stats_final']).reset_index(drop=True)
            data_mean_subject['sub_origine'] = data['subject'].unique()[0]

            # Keep te rois on the good order 
            data_mean_subject['rois'] = pd.Categorical(data_mean_subject['rois'], categories=rois_list, ordered=True)
            data_mean_subject = data_mean_subject.sort_values(by='rois')
            data_mean_all = pd.concat([data_mean_all, data_mean_subject], ignore_index=True)

        # Rename columns and set 'group' as subject
        data_all = data_all.rename(columns={'subject': 'sub-origine'})
        data_all['subject'] = ['group'] * len(data_all)
        data_mean_all['subject'] = ['group'] * len(data_mean_all)

        # Export aggregated data to TSV
        tsv_all_dir = '{}/{}/derivatives/pp_data/group/{}/prf/tsv'.format(main_dir, project_dir, format_)
        os.makedirs(tsv_all_dir, exist_ok=True)
        data_all.to_csv('{}/group_css-prf_{}.tsv'.format(tsv_all_dir, data_type), sep="\t", na_rep='NaN', index=False)
        data_mean_all.to_csv('{}/group_css-prf_{}.tsv'.format(tsv_all_dir, data_type), sep="\t", na_rep='NaN', index=False)



# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))
  
    

    
    






