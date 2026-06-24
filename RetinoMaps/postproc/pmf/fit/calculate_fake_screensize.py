"""
-----------------------------------------------------------------------------------------
calculate_fake_screensize.py
-----------------------------------------------------------------------------------------
Goal of the script:
Calculate fake screen size to accomodate bigger design matrix in pMF analysis
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
python calculate_fake_screensize.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""

import json
import os
import math
import yaml
import sys

# --- Config ---
main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]

sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

# --- Load settings ---
base_dir               = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path  = os.path.join(base_dir, project_dir, 'settings.yml')
analysis_settings_path = os.path.join(base_dir, project_dir, 'pmf-analysis.yml')
analysis_info          = load_settings([general_settings_path, analysis_settings_path])[0]

task     = analysis_info["analysis_task_names"][0]
subjects = ["sub-01"]

# --- Visual angle settings ---
original_dva = 20   # degrees: original stimulus window
fake_dva     = 30   # degrees: new (faked) stimulus window

# --- YML output path (on scratch) ---
analysis_settings_path = os.path.join(main_dir, project_dir, 'pmf-analysis.yml')

# --- Find first valid JSON to read screen info ---
# (screen settings are the same across subjects, so we just need one)
screen_info = None

for subject in subjects:
    for run in range(0, 2):
        session   = "ses-01" if subject == "sub-01" else "ses-02"
        json_path = f"{main_dir}/{project_dir}/{subject}/{session}/func/{subject}_{session}_task-{task}_run-0{run+1}_events.json"

        if not os.path.exists(json_path):
            continue

        with open(json_path, "r") as f:
            data = json.load(f)

        stim = data['StimulusPresentation']
        screen_info = {
            'distance_m' : stim['ScreenDistance'],   # already in meters
            'size_m'     : stim['ScreenSize'],        # [x, y] in meters
            'resolution' : stim['ScreenResolution'],  # [x, y] in pixels
        }
        print(f"Screen info read from: {json_path}")
        break
    if screen_info:
        break

if screen_info is None:
    raise RuntimeError("No valid JSON file found — check paths.")

# --- Unpack (convert to cm for geometry) ---
distance_cm      = screen_info['distance_m'] * 100
size_x_cm        = screen_info['size_m'][0]  * 100
size_y_cm        = screen_info['size_m'][1]  * 100
res_x, res_y     = screen_info['resolution']

# --- Pixels per degree from real screen geometry ---
# Full screen subtends: dva = 2 * atan(half_screen_cm / distance_cm)
screen_dva_x = 2 * math.degrees(math.atan((size_x_cm / 2) / distance_cm))
screen_dva_y = 2 * math.degrees(math.atan((size_y_cm / 2) / distance_cm))
pix_per_deg  = res_x / screen_dva_x   # px per degree of visual angle

print(f"\nOriginal setup:")
print(f"  Resolution        : {res_x} x {res_y} px")
print(f"  Screen size       : {size_x_cm:.4f} x {size_y_cm:.4f} cm")
print(f"  Distance          : {distance_cm:.4f} cm")
print(f"  Screen size       : {screen_dva_x:.2f} x {screen_dva_y:.2f} dva (full screen)")
print(f"  Stimulus window   : {original_dva} dva")
print(f"  Pixels/degree     : {pix_per_deg:.4f}")

# --- Fake resolution: how many pixels span fake_dva at the same px/deg? ---
fake_res_x = round(pix_per_deg * fake_dva)
fake_res_y = round(pix_per_deg * fake_dva)  # square window

# --- Fake screen size in cm: physical width that subtends fake_dva at same distance ---
# size_cm = 2 * distance_cm * tan(fake_dva/2)
fake_size_cm_x = 2 * distance_cm * math.tan(math.radians(fake_dva / 2))
fake_size_cm_y = fake_size_cm_x  # square window

print(f"\nFake setup (scaled to {fake_dva} dva):")
print(f"  Fake resolution   : {fake_res_x} x {fake_res_y} px")
print(f"  Fake size         : {fake_size_cm_x:.4f} x {fake_size_cm_y:.4f} cm")

# --- Load existing yml (if present) and update ---
if os.path.exists(analysis_settings_path):
    with open(analysis_settings_path, "r") as f:
        settings = yaml.safe_load(f) or {}
else:
    settings = {}

settings['fake_screen_size_pix'] = [fake_res_x, fake_res_y]
settings['fake_screen_size_cm']  = [round(fake_size_cm_x, 4), round(fake_size_cm_y, 4)]

with open(analysis_settings_path, "w") as f:
    yaml.dump(settings, f, default_flow_style=False)

print(f"\nSaved to: {analysis_settings_path}")