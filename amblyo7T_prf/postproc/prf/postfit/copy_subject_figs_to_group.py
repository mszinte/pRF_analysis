"""
-----------------------------------------------------------------------------------------
copy_subject_figs_to_group.py
-----------------------------------------------------------------------------------------
Goal of the script:
Copy individual subject figures (correlation, active vertices, ecc-size-pcm) to group
folders, organizing them by figure type and task in subdirectories.

For each group (group-patient, group-control):
- Loads list of subjects from settings
- Copies correlation figures to group-patient/fsnative/prf/figures-corr/task-name/
- Copies active vertices figures to group-patient/fsnative/prf/figures-active-vert/task-name/
- Copies ecc-size-pcm figures to group-patient/fsnative/prf/figures-ecc-size-pcm/task-name/
- Renames files to include subject name: subject_figname.pdf
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /scratch/mszinte/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: group label (e.g. group-patient or group-control)
-----------------------------------------------------------------------------------------
Output(s):
Figures copied to group folders organized by figure type and task
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python copy_subject_figs_to_group.py [main directory] [project name] [group label]

Examples:
python copy_subject_figs_to_group.py /scratch/mszinte/data amblyo7T_prf group-patient
python copy_subject_figs_to_group.py /scratch/mszinte/data amblyo7T_prf group-control
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import shutil
import re
import subprocess
import tempfile

# Personal import
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group_label = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

# Get subjects for this group
if 'patient' in group_label:
    subjects_to_group = analysis_info['group_patient']
elif 'control' in group_label:
    subjects_to_group = analysis_info['group_control']
else:
    raise ValueError(f"Unknown group label: {group_label}")

formats = analysis_info['formats']

# Figure types and their subdirectories
figure_types = {
    '_corr': 'figures-corr',
    '_active-vert': 'figures-active-vert',
    '_ecc-size-pcm': 'figures-ecc-size-pcm'
}

# Process each format
for format_ in formats:
    
    # Copy figures from each subject
    for subject in subjects_to_group:
        print(f'Processing figures for {subject} - {format_}')
        
        subject_fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
            main_dir, project_dir, subject, format_)
        
        if not os.path.isdir(subject_fig_dir):
            print(f"[SKIP] Figure directory not found: {subject_fig_dir}")
            continue
        
        # Copy each figure type
        for fig_type, fig_subdir in figure_types.items():
            
            # Find figures matching this type
            for fn in os.listdir(subject_fig_dir):
                
                if fig_type in fn:
                    # Skip eye-specific figures (RightEye/LeftEye)
                    if 'RightEye' in fn or 'LeftEye' in fn:
                        continue
                    
                    # Extract task name from filename (e.g., BarsBars from task-BarsBars_fmriprep...)
                    match = re.search(r'task-([A-Za-z]+)_', fn)
                    if match:
                        task_name = match.group(1)
                    else:
                        task_name = 'no-task'
                    
                    # Create task subfolder within figure type folder
                    group_fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/{}/{}'.format(
                        main_dir, project_dir, group_label, format_, fig_subdir, task_name)
                    os.makedirs(group_fig_dir, exist_ok=True)
                    
                    src = os.path.join(subject_fig_dir, fn)
                    # Rename with subject prefix
                    new_fn = f'{subject}_{fn}'
                    dst = os.path.join(group_fig_dir, new_fn)
                    
                    # Use rsync to copy with temp name, then rename
                    with tempfile.NamedTemporaryFile(dir=group_fig_dir, delete=False) as tmp:
                        tmp_path = tmp.name
                    
                    try:
                        # rsync with -a flag (archive mode: preserves permissions, timestamps, etc)
                        subprocess.run(['rsync', '-a', src, tmp_path], 
                                     check=True, 
                                     stdout=subprocess.DEVNULL,
                                     stderr=subprocess.DEVNULL)
                        # Rename to final name
                        os.rename(tmp_path, dst)
                        print(f'  Copying: {fn}')
                    except subprocess.CalledProcessError as e:
                        print(f'  [ERROR] rsync failed for {fn}: {e}')
                        if os.path.exists(tmp_path):
                            os.remove(tmp_path)
                    except Exception as e:
                        print(f'  [ERROR] Failed to copy {fn}: {e}')
                        if os.path.exists(tmp_path):
                            os.remove(tmp_path)

print(f'Done copying figures to {group_label}')