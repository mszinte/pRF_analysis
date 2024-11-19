def extract_predictions_r2 (labels,estimate,source_data):
    """
    extract prediction and correponding r2 from Nilearn run_glm outputs
    
    Parameters
    ----------
    labels : labels outptuts from Nilearn run_glm
    estimate : estimate outptuts from Nilearn run_glm
    source_data : data uses for run_glm
    
    Returns
    -------
    predictions : the predictions from the glm 
    rsquare : the rsqare of the predictions 
    
    """
    
    import numpy as np
    cl_names = np.unique(labels)

    n_cl = cl_names.size
    list_cl = []

    idx_in_cl = np.zeros_like(labels, dtype=int)

    for i in range(n_cl):
        list_idx = np.argwhere(labels==cl_names[i]).flatten()
        list_cl += [list_idx]
        for i_cnt, idx in enumerate(list_idx):
            idx_in_cl[idx] = i_cnt

    # acces to predictions and r2 for left hemisphere 
    predictions = np.zeros(source_data.shape, dtype=float)
    rsquare = np.zeros_like(labels, dtype=float)
    N = rsquare.size
    for idx,labels_ in enumerate(labels) :    
        label_mask = labels[idx]
        predictions[:,idx] = estimate[label_mask].predicted[:,idx_in_cl[idx]]
        rsquare[idx] = estimate[label_mask].r_square[idx_in_cl[idx]]
        
    return predictions, rsquare

def eventsMatrix(design_file, task, tr):
    """
    Returns the events matrix for the GLM. Works for the Sac/Pur Localisers and Sac/Pur Visual/Endogenous Localisers
    Parameters
    ----------
    design_file         : path to the tsv/csv file containing the events
    task                : task name (SacLoc / PurLoc / SacVELoc / PurVELoc)
    tr                  : time repetition in seconds (e.g. 1.2)
    Returns
    -------
    new_events_glm - pandas DataFrame containing the events for the GLM 
    """
    import pandas as pd
    import numpy as np
    import ipdb
    deb = ipdb.set_trace
    tr_dur = tr
    events = pd.read_table(design_file)


    if 'pMF' in task:
        events_glm = events[['onset','duration','trial_type']].copy(deep=True)
        events_glm.replace({'trial_type': {3: 'Fix', 1: 'PurSac', 2: 'PurSac'}},inplace=True)
        events_glm.duration = tr_dur
        events_glm.onset = 0
        
        events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())
        
        new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
        
        for idx, group in enumerate(events_glm_groups):
            onset = np.round(group[1]['onset'][group[1].index[0]],2)
            dur = np.round(sum(group[1]['duration']),2)
            ttype = group[1]['trial_type'][group[1].index[0]]
            
            new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'trial_type'])
            new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)
            
    elif 'pRF' in task:
        events_glm = events[['onset','duration','bar_direction']].copy(deep=True)
        events_glm = events[['onset', 'duration', 'bar_direction']].replace({'bar_direction': {1: 'vision', 3: 'vision', 5: 'vision', 7: 'vision', 9: 'Fix'}}).rename(columns={'bar_direction': 'trial_type'})
        
        events_glm.onset = 0
        
        events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())
        
        new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
        
        for idx, group in enumerate(events_glm_groups):
            onset = np.round(group[1]['onset'][group[1].index[0]],2)
            dur = np.round(sum(group[1]['duration']),2)
            ttype = group[1]['trial_type'][group[1].index[0]]
        
            new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'trial_type'])
            new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)

    else:
        if 'VE' not in task: # Sac/Pur Loc
            events_glm = events[['onset','duration','trial_type']].copy(deep=True)
            events_glm.replace({'trial_type': {3: 'Fix', 1: 'Sac', 2: 'Pur'}},inplace=True)
            events_glm.duration = tr_dur
            events_glm.onset = 0
            events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())

            new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
            for idx, group in enumerate(events_glm_groups):
                onset = np.round(group[1]['onset'][group[1].index[0]],2)
                dur = np.round(sum(group[1]['duration']),2)
                ttype = group[1]['trial_type'][group[1].index[0]]

                new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'trial_type'])
                new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)


        elif 'SacVE' in task: # Sac/Pur Loc: # Visual-Endogenous Sac/Pur Loc 
            ###### OLD way ####
            # events_glm = events[['onset','duration','trial_type', 'eyemov_vis_end']].copy(deep=True)
            # events_glm.replace({'trial_type': {3: 'Fix', 1: 'Sac', 2: 'Pur'}},inplace=True)
            # events_glm.replace({'eyemov_vis_end': {3: 'Fix', 1: 'SacExo', 2: 'SacEndo'}},inplace=True)
            # events_glm.duration = tr_dur
            # events_glm.onset = 0
            # events_glm_groups = events_glm.groupby((events_glm.eyemov_vis_end!=events_glm.eyemov_vis_end.shift()).cumsum())

            # new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
            
            # for idx, group in enumerate(events_glm_groups):
            #     onset = np.round(group[1]['onset'][group[1].index[0]],2)
            #     dur = np.round(sum(group[1]['duration']),2)
            #     ttype = group[1]['trial_type'][group[1].index[0]]
            #     vis_end = group[1]['eyemov_vis_end'][group[1].index[0]]

            #     new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'vis_end'])
            #     new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)
            
            ###### New way ####
            events_glm = events[['onset','duration','eyemov_vis_end']].copy(deep=True)
            events_glm = events[['onset', 'duration', 'eyemov_vis_end']].replace({'eyemov_vis_end': {3: 'Fix', 1: 'SacExo', 2: 'SacEndo'}}).rename(columns={'eyemov_vis_end': 'trial_type'})
            
            events_glm.onset = 0
            
            events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())
            
            new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
            
            for idx, group in enumerate(events_glm_groups):
                onset = np.round(group[1]['onset'][group[1].index[0]],2)
                dur = np.round(sum(group[1]['duration']),2)
                ttype = group[1]['trial_type'][group[1].index[0]]
            
                new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'trial_type'])
                new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)
                        

        elif 'PurVE' in task: # Sac/Pur Loc: # Visual-Endogenous Sac/Pur Loc 
            # events_glm = events[['onset','duration','trial_type', 'eyemov_vis_end']].copy(deep=True)
            # events_glm.replace({'trial_type': {3: 'Fix', 1: 'Sac', 2: 'Pur'}},inplace=True)
            # events_glm.replace({'eyemov_vis_end': {3: 'Fix', 1: 'PurExo', 2: 'PurEndo'}},inplace=True)
            # events_glm.duration = tr_dur
            # events_glm.onset = 0
            # events_glm_groups = events_glm.groupby((events_glm.eyemov_vis_end!=events_glm.eyemov_vis_end.shift()).cumsum())

            # new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type', 'vis_end'])
            # for idx, group in enumerate(events_glm_groups):
            #     onset = np.round(group[1]['onset'][group[1].index[0]],2)
            #     dur = np.round(sum(group[1]['duration']),2)
            #     ttype = group[1]['trial_type'][group[1].index[0]]
            #     vis_end = group[1]['eyemov_vis_end'][group[1].index[0]]

            #     new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'vis_end'])
            #     new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)
            events_glm = events[['onset','duration','eyemov_vis_end']].copy(deep=True)
            events_glm = events[['onset', 'duration', 'eyemov_vis_end']].replace({'eyemov_vis_end': {3: 'Fix', 1: 'PurExo', 2: 'PurEndo'}}).rename(columns={'eyemov_vis_end': 'trial_type'})
            
            events_glm.onset = 0
            
            events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())
            
            new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
            
            for idx, group in enumerate(events_glm_groups):
                onset = np.round(group[1]['onset'][group[1].index[0]],2)
                dur = np.round(sum(group[1]['duration']),2)
                ttype = group[1]['trial_type'][group[1].index[0]]
            
                new_row = pd.Series([onset, dur, ttype], index=['onset', 'duration', 'trial_type'])
                new_events_glm = pd.concat([new_events_glm, new_row.to_frame().T], ignore_index=True)

    for idx in new_events_glm.index:
        if idx==0:
            new_events_glm.at[idx, 'onset'] = 0
        else:
            new_events_glm.at[idx, 'onset'] = np.round(new_events_glm.at[idx-1, 'onset'] + new_events_glm.at[idx-1, 'duration'],2)

    return new_events_glm


def make_design(events, tr, n_times, hrf=None):
    """generate either a blip design or one convolved with an hrf

    Args:
        events ([type]): [description]
        tr ([type]): [description]
        n_times ([type]): [description]
        hrf ([type], optional): Defaults to None. [description]

    Returns:
        [type]: [description]
    """
    import numpy as np

    # loop over conditions
    conditions = np.unique(events.trial_type)
    n_conditions = len(set(events['trial_type'].values))

    dm = np.zeros((n_times, n_conditions))

    if hrf is None:
        for i, cond in enumerate(conditions):

            # onset times for qth condition in run p
            otimes = np.array(
                (events[events['trial_type'] == cond]['onset'].values//tr)).astype(int)
            yvals = np.zeros((n_times))
            for r in otimes:
                yvals[r] = 1
            dm[:, i] = yvals

    else:
        # calc
        all_times = np.linspace(0, tr*(n_times-1), n_times)
        hrf_times = np.linspace(0, tr*(len(hrf)-1), len(hrf))

        for i, cond in enumerate(conditions):
            # onset times for qth condition in run p
            otimes = events[events['trial_type'] == cond]['onset'].values

            # intialize
            yvals = np.zeros((n_times))

            # loop over onset times
            for r in otimes:
                # interpolate to find values at the data sampling time points
                f = pchip(r + hrf_times, hrf, extrapolate=False)(all_times)
                yvals = yvals + np.nan_to_num(f)

            # record
            dm[:, i] = yvals

    return dm