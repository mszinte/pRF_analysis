"""
-----------------------------------------------------------------------------------------
get_screen_info_mat.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate motion design matrix and videos of saccade direction in retinal space.
Uses saccade analysis from RetinoMaps. 
Loops over all runs and creates concatenated downsampled MDM.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Updated event json files
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python get_screen_info_mat.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""

import scipy.io
import json
import os
import sys

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]

sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

# ── Load settings + events 
base_dir      = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, 'settings.yml')
analysis_settings_path = os.path.join(base_dir, project_dir, 'pmf-analysis.yml')
analysis_info = load_settings([general_settings_path, analysis_settings_path])[0]


task = analysis_info["analysis_task_names"][0]
subjects = analysis_info["subjects"]

for subject in subjects:
    for run in range(0,2): 
        session = "ses-01" if subject == "sub-01" else "ses-02"

        # --- Paths ---
        source_data_path = f"{main_dir}/{project_dir}/sourcedata/{subject}/{session}/add"
        matfile = f"{source_data_path}/{subject}_{session}_task-{task}_run-0{run+1}_matFile.mat"
        json_path = f"{main_dir}/{project_dir}/{subject}/{session}/func/{subject}_{session}_task-{task}_run-0{run+1}_events.json"

        # --- Load .mat file ---
        if not os.path.exists(matfile):
            print(f"[SKIP] Mat file not found: {matfile}")
            continue

        print(f"\n{'='*60}")
        print(f"Processing {subject} | {task} | {session} | run-0{run+1}")

        mat = scipy.io.loadmat(matfile)
        config = mat['config']
        cfg0 = config[0][0]
        scr = cfg0['scr'][0][0]

        # --- Extract screen info (stored in cm, convert to m) ---
        screen_res_x  = int(scr['scr_sizeX'][0][0])
        screen_res_y  = int(scr['scr_sizeY'][0][0])
        screen_size_x = float(scr['disp_sizeX'][0][0]) / 1000  
        screen_size_y = float(scr['disp_sizeY'][0][0]) / 1000   
        distance_m    = float(scr['dist'][0][0]) / 100         

        print(f"  Screen resolution : {screen_res_x} x {screen_res_y} px")
        print(f"  Screen size       : {screen_size_x:.4f} x {screen_size_y:.4f} m")
        print(f"  Viewing distance  : {distance_m:.4f} m")

        # --- Load and update JSON ---
        if not os.path.exists(json_path):
            print(f"[SKIP] JSON file not found: {json_path}")
            continue

        with open(json_path, "r") as f:
            data = json.load(f)

        data['StimulusPresentation']['ScreenDistance']   = distance_m
        data['StimulusPresentation']['ScreenSize']       = [screen_size_x, screen_size_y]
        data['StimulusPresentation']['ScreenResolution'] = [screen_res_x, screen_res_y]

        with open(json_path, "w") as f:
            json.dump(data, f, indent=4)

        print(f"  JSON updated: {json_path}")

print(f"\nSetting permissions in {main_dir}/{project_dir}")
os.system(f"chmod -Rf 771 {main_dir}/{project_dir}")
os.system(f"chgrp -Rf {group} {main_dir}/{project_dir}")
print("Done.")