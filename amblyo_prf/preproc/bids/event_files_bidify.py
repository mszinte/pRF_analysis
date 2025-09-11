"""
-----------------------------------------------------------------------------------------
event_files_bidify.py
-----------------------------------------------------------------------------------------
Goal of the script:
* Transform the original event files in a tsv bids format
-----------------------------------------------------------------------------------------
Input(s) in that order:
* path to the project data folder (stereo_prf) - inside that folder, you should find sourcedata, derivatives, ...
* -o : overwrite optional parameter: -o for overwriting, otherwise skip
-----------------------------------------------------------------------------------------
Output(s):
the BIDS tsv event file correctly placed in the output structure in dcm2niix folder
-----------------------------------------------------------------------------------------
To run:
Find where event_files_bidify.py is located (say in ~/XX/), and where your project files are located (say ~/path2stereo_prf/) and run:
>> python ~/XX/event_files_bidify.py ~/path2stereo_prf/
Example:
From mesocentre:
>> python  ~/projects/stereo_prf/analysis_code/preproc/bids/event_files_bidify.py '/scratch/mszinte/data/stereo_prf/' 
From invibe nohost:
>> python  ~/disks/meso_H/projects/stereo_prf/analysis_code/preproc/bids/event_files_bidify.py '~/disks/meso_S/data/stereo_prf/' -o
-----------------------------------------------------------------------------------------
Written by Adrien Chopin (adrien.chopin@gmail.com) - 2022
-----------------------------------------------------------------------------------------
"""

import pandas as pd
import os
import sys

project_dir = os.path.expanduser(sys.argv[1])

# define overwriting (or skipping) behavior
overwrite = 0            # default skip
if len(sys.argv)>2:      # if more than two arguments, then it must be the overwrite parameter
    if sys.argv[2]=="-o":
        overwrite = 1    # it is!

# define data directories, roots for source and destination
rootpath = os.path.join(project_dir,'sourcedata') 
sourceroot = os.path.join(rootpath,'Big_data_STAM','prf_event_files')
destroot = os.path.join(rootpath,'dcm2niix')

bar_dic = {" LL_UR ": "1"," Right ":"2", " LR_UL ":"3", " Down ":"4", " Up ":"5", " UL_LR ":"6", " Left ":"7", " UR_LL ":"8","LL_UR": "1","Right":"2", "LR_UL":"3", "Down":"4", "Up":"5", "UL_LR":"6", "Left":"7", "UR_LL":"8"}
event_dic = {' fixation dimming ':'1', ' hit ':'2', ' miss ':'3', ' false_alarm ':'4','fixation dimming':'1', 'hit':'2', 'miss':'3', 'false_alarm':'4'}

for file in os.listdir(sourceroot):
    if file.endswith('.tsv'):
        filepath = os.path.join(sourceroot,file)
        print('Working with '+file)
        
        # define destination file       
        dest_path = os.path.join(destroot,file[0:6],file[7:13],'func',file[0:28]+file[33:44])
       # dest_path2 = os.path.join(project_dir,file[0:6],file[7:13],'func',file[0:28]+file[33:44])

        # overwriting or skipping
        if os.path.exists(dest_path) & (overwrite==0):
            print('File already exists - skip '+file)
        else:
            if os.path.exists(dest_path) & (overwrite==1):
                print('File already exists - overwrite '+dest_path) 

            # import file renaming cols
            data=pd.read_csv(filepath,sep='\t',names=['onset', 'bar_direction_str','event_type_str'],header=0) 
            data['trial_number'] = data.index+1   #create a trial number col with index

            # transform categorical col to a numbered col
            # Do that by remaping the values of the bar_direction_str col using a dictionary
            data['bar_direction'] = data['bar_direction_str'].map(bar_dic)
            data.drop('bar_direction_str', inplace=True, axis=1) # and remove old col

            # Remap the values of the event_type_str col
            data['event_type'] = pd.to_numeric(data['event_type_str'].map(event_dic))
            data.drop('event_type_str', inplace=True, axis=1) # and remove old col

            # false alarms events are instants, max 100 ms (we correct the others to 0.8 sec below)
            data['duration'] = 0.1 
            # for all events, except false alarms, the event actually happened 0.8 sec before
            data.loc[data['event_type']<4,'onset'] = data.loc[data['event_type']<4,'onset'] - 0.8
            data.loc[data['event_type']<4,'duration'] = 0.8
           
            #change order of columns
            cols = ['onset','duration','trial_number','bar_direction','event_type']                    
            data = data[cols]
            
            # save tsv file
            print('Saving to '+dest_path)
            data.to_csv(dest_path, index = False, sep='\t')
            #print('Saving to '+dest_path2)
            #data.to_csv(dest_path2, index = False, sep='\t')
print('Event file conversion done')


