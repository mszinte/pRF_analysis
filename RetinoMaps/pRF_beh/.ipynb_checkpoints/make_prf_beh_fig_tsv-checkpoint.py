"""
-----------------------------------------------------------------------------------------
make_prf_beh_fig_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
make tsv for pRF behaviour figures 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
pRF behaviour analysis tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/pRF_beh/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/pRF_beh/
python make_prf_beh_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_prf_beh_fig_tsv.py /scratch/mszinte/data RetinoMaps group 327
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
import json
import numpy as np
import pandas as pd

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
n_runs_prf = analysis_info['n_runs_prf']
subject_excluded = analysis_info['outlier_prf_beh']

if 'group' not in subject:
    print('{} is processing...'.format(subject))
    # Exception for sub-01
    if subject == 'sub-01': session = 'ses-02'
    else: session = 'ses-01'
    for n_run in range(n_runs_prf):
        prf_events_dir = '{}/{}/{}/{}/func'.format(main_dir, project_dir, subject, session)
        prf_events_fn = '{}_{}_task-pRF_run-0{}_events.tsv'.format(subject, session ,n_run+1)
        prf_events_raw = pd.read_table('{}/{}'.format(prf_events_dir, prf_events_fn))
        
        # Find TRs with bar stimulation
        prf_events = prf_events_raw.loc[prf_events_raw['bar_direction'] != 9].copy()
            
        # concat runs
        if n_run==0: 
            perf_runs_indiv_concat = prf_events.response_val.values
            stair_runs_indiv_concat = prf_events.stim_stair_val.values
        else: 
            perf_runs_indiv_concat = np.concatenate([perf_runs_indiv_concat, prf_events.response_val.values])
            stair_runs_indiv_concat = np.concatenate([stair_runs_indiv_concat, prf_events.stim_stair_val.values])
            
    # Compute cumulative mean of median and staircase
    perf_runs_indiv = np.zeros(stair_runs_indiv_concat.shape[0])
    stair_runs_indiv = np.zeros(stair_runs_indiv_concat.shape[0])
    for t in range(stair_runs_indiv_concat.shape[0]):
                    perf_runs_indiv[t] = np.nanmean(perf_runs_indiv_concat[0:t])
                    stair_runs_indiv[t] = np.nanmean(stair_runs_indiv_concat[0:t])
                    
    # Make DF 
    # Trials DF
    prf_beh_trials_df_indiv = pd.DataFrame({
        'subject':[subject]*perf_runs_indiv.shape[0], 
        'perf': perf_runs_indiv, 
        'stair': stair_runs_indiv})
    
    # median DF
    prf_beh_median_df_indiv = pd.DataFrame()
    prf_beh_median_df_indiv['subject'] = [subject]
    prf_beh_median_df_indiv['perf_median'] = [np.nanmedian(perf_runs_indiv)]
    prf_beh_median_df_indiv['perf_ci_low'] = [np.nanquantile(perf_runs_indiv, 0.025)]
    prf_beh_median_df_indiv['perf_ci_up'] = [np.nanquantile(perf_runs_indiv, 0.975)]
    prf_beh_median_df_indiv['stair_median'] = [np.nanmedian(stair_runs_indiv)]
    prf_beh_median_df_indiv['stair_ci_low'] = [np.nanquantile(stair_runs_indiv, 0.025)]
    prf_beh_median_df_indiv['stair_ci_up'] = [np.nanquantile(stair_runs_indiv, 0.975)] 
    
    # Export DF as TSV 
    prf_beh_tsv_dir = '{}/{}/derivatives/pp_data/{}/pRF_beh/tsv'.format(main_dir, project_dir, subject)
    os.makedirs(prf_beh_tsv_dir, exist_ok=True)
    prf_beh_trials_tsv_fn = '{}/{}_pRF_beh_trials.tsv'.format(prf_beh_tsv_dir, subject)
    prf_beh_median_tsv_fn = '{}/{}_pRF_beh_median.tsv'.format(prf_beh_tsv_dir, subject)

    print('export {} and {}'.format(prf_beh_trials_tsv_fn, prf_beh_median_tsv_fn))
    prf_beh_trials_df_indiv.to_csv(prf_beh_trials_tsv_fn, sep="\t", na_rep='NaN', index=False)
    prf_beh_median_df_indiv.to_csv(prf_beh_median_tsv_fn, sep="\t", na_rep='NaN', index=False)

else :
    # group processing
    print('{} is processing...'.format(subject))
    for n_subject, subject in enumerate(subjects):
        prf_beh_tsv_dir = '{}/{}/derivatives/pp_data/{}/pRF_beh/tsv'.format(main_dir, project_dir, subject)
        prf_beh_trials_tsv_fn = '{}/{}_pRF_beh_trials.tsv'.format(prf_beh_tsv_dir, subject)
        prf_beh_median_tsv_fn = '{}/{}_pRF_beh_median.tsv'.format(prf_beh_tsv_dir, subject)
        
        prf_beh_trials_df_indiv = pd.read_table(prf_beh_trials_tsv_fn)
        prf_beh_median_df_indiv = pd.read_table(prf_beh_median_tsv_fn)
        
        # Concat DF from subjects
        if n_subject == 0: 
            prf_beh_trials_df_group = prf_beh_trials_df_indiv.copy()
            prf_beh_median_df_group = prf_beh_median_df_indiv.copy()
        else: 
            prf_beh_trials_df_group = pd.concat([prf_beh_trials_df_group, prf_beh_trials_df_indiv])
            prf_beh_median_df_group = pd.concat([prf_beh_median_df_group, prf_beh_median_df_indiv])
            
    # Export DF as TSV 
    prf_beh_tsv_dir = '{}/{}/derivatives/pp_data/group/pRF_beh/tsv'.format(main_dir, project_dir)
    os.makedirs(prf_beh_tsv_dir, exist_ok=True)
    prf_beh_trials_group_tsv_fn = '{}/group_pRF_beh_trials.tsv'.format(prf_beh_tsv_dir)
    prf_beh_median_group_tsv_fn = '{}/group_pRF_beh_median.tsv'.format(prf_beh_tsv_dir)
    
    print('export {} and {}'.format(prf_beh_trials_group_tsv_fn, prf_beh_median_group_tsv_fn))
    prf_beh_trials_df_group.to_csv(prf_beh_trials_group_tsv_fn, sep="\t", na_rep='NaN', index=False)
    prf_beh_median_df_group.to_csv(prf_beh_median_group_tsv_fn, sep="\t", na_rep='NaN', index=False)
    
    # Make inter subject median DF without outliers 
    
    # Trials DF 
    # Exclude outliers
    prf_beh_trials_df_group_filtered = prf_beh_trials_df_group[~prf_beh_trials_df_group['subject'].isin(subject_excluded)]
    
    # compute inter subject median and CI 
    perf_pivot = prf_beh_trials_df_group_filtered.pivot(columns='subject', values='perf')
    stair_pivot = prf_beh_trials_df_group_filtered.pivot(columns='subject', values='stair')
    
    
    median_perf = perf_pivot.median(axis=1)
    ci_low_perf = perf_pivot.quantile(0.025, axis=1)
    ci_up_perf = perf_pivot.quantile(0.975, axis=1)
    
    median_stair = stair_pivot.median(axis=1)
    ci_low_stair = stair_pivot.quantile(0.025, axis=1)
    ci_up_stair = stair_pivot.quantile(0.975, axis=1)
    
    # Creat DF for inter subject median and CI
    prf_beh_trials_df_median_group = pd.DataFrame({
        'subject': ['group'] * len(median_perf),
        'median_perf': median_perf,
        'ci_low_perf': ci_low_perf,
        'ci_up_perf': ci_up_perf,
        'median_stair': median_stair,
        'ci_low_stair': ci_low_stair,
        'ci_up_stair': ci_up_stair,
    })
    
    # Add outlier to DF
    prf_beh_trials_df_group_excluded = prf_beh_trials_df_group[prf_beh_trials_df_group['subject'].isin(subject_excluded)]
    prf_beh_trials_df_group_excluded = prf_beh_trials_df_group_excluded.rename(columns={'perf': 'median_perf', 'stair': 'median_stair'})
 
    # Concat DF
    prf_beh_trials_df_median_group = pd.concat([prf_beh_trials_df_median_group, prf_beh_trials_df_group_excluded], ignore_index=True)
    
    # Median DF
    # Exclude outliers
    prf_beh_median_df_group_filtered = prf_beh_median_df_group[~prf_beh_median_df_group['subject'].isin(subject_excluded)]
    prf_beh_median_df_group_excluded = prf_beh_median_df_group[prf_beh_median_df_group['subject'].isin(['sub-07', 'sub-12'])]
    
    # Cpmpute median
    median_values = prf_beh_median_df_group_filtered.drop(columns=['subject']).median()
    
    # concat DF 
    prf_beh_df_median_subject = pd.DataFrame(median_values).T
    prf_beh_df_median_subject['subject'] = 'group'
    prf_beh_df_median_subject = pd.concat([prf_beh_df_median_subject, prf_beh_median_df_group_excluded])
    
    # Export DF as TSV 
    prf_beh_trials_group_median_tsv_fn = '{}/group_pRF_beh_trials_group_median.tsv'.format(prf_beh_tsv_dir)
    prf_beh_median_group_median_tsv_fn = '{}/group_pRF_beh_median_group_median.tsv'.format(prf_beh_tsv_dir)
    
    print('export {} and {}'.format(prf_beh_trials_group_median_tsv_fn, prf_beh_median_group_median_tsv_fn))
    prf_beh_trials_df_median_group.to_csv(prf_beh_trials_group_median_tsv_fn, sep="\t", na_rep='NaN', index=False)
    prf_beh_df_median_subject.to_csv(prf_beh_median_group_median_tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))  

        
       