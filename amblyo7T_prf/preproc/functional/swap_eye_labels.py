"""
-----------------------------------------------------------------------------------------
swap_eye_labels.py
-----------------------------------------------------------------------------------------
Goal of the script:
Swap LeftEye <-> RightEye labels in fMRIPrep output
Preserves entire filename structure, backs up original func/ to func_incorrect_label/
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: session name
sys.argv[5]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# Renamed fMRIPrep output files with corrected eye labels
# Original files backed up to func_incorrect_label/
-----------------------------------------------------------------------------------------
To run:
python swap_eye_labels.py [main directory] [project name] [subject] [session] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/functional/
python swap_eye_labels.py /scratch/mszinte/data amblyo7T_prf sub-10 ses-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
with Claude (claude.ai)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import shutil
import datetime
from pathlib import Path

# Time
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
session = sys.argv[4]
group = sys.argv[5]

# Construct paths
sub_dir = Path(main_dir) / project_dir / "derivatives" / "fmriprep" / "fmriprep" / subject
func_dir = sub_dir / session / "func"
backup_dir = sub_dir / session / "func_incorrect_label"

# Verify func directory exists
if not func_dir.exists():
    print(f"Error: {func_dir} does not exist")
    sys.exit(1)

print(f"Processing: {subject} / {session}")
print(f"Source directory: {func_dir}")
print(f"Backup directory: {backup_dir}")
print()

# Count files
files = list(func_dir.glob("*"))
left_eye_files = [f for f in files if f.is_file() and "LeftEye" in f.name]
right_eye_files = [f for f in files if f.is_file() and "RightEye" in f.name]

print(f"Found {len(left_eye_files)} files with LeftEye")
print(f"Found {len(right_eye_files)} files with RightEye")
print()

if len(left_eye_files) == 0 and len(right_eye_files) == 0:
    print("Warning: No files with LeftEye or RightEye found. Exiting.")
    sys.exit(0)

# Step 1: Backup original func directory
if backup_dir.exists():
    print(f"Warning: {backup_dir} already exists. Skipping backup.")
else:
    print("Step 1: Creating backup...")
    func_dir.rename(backup_dir)
    print(f"✓ Backed up to: {backup_dir}")

# Step 2: Create new func directory and swap files
print()
print("Step 2: Swapping eye labels and copying files...")
func_dir.mkdir(parents=True, exist_ok=True)

for old_file in backup_dir.glob("*"):
    if not old_file.is_file():
        continue
    
    filename = old_file.name
    
    if "LeftEye" in filename:
        # LeftEye -> RightEye
        newname = filename.replace("LeftEye", "RightEye")
        print(f"  {filename} -> {newname}")
        
    elif "RightEye" in filename:
        # RightEye -> LeftEye
        newname = filename.replace("RightEye", "LeftEye")
        print(f"  {filename} -> {newname}")
        
    else:
        # No eye label, copy as-is
        newname = filename
    
    # Use rsync for faster copying with progress
    src = str(old_file)
    dst = str(func_dir / newname)
    os.system(f"rsync -av --progress {src} {dst}")

# Step 3: Handle figures directory
figures_dir = sub_dir / "figures"
figures_backup_dir = sub_dir / "figures_incorrect_label"

if figures_dir.exists():
    print()
    print("Step 3: Processing figures directory...")
    
    fig_files = list(figures_dir.glob("*"))
    left_eye_figs = [f for f in fig_files if f.is_file() and "LeftEye" in f.name]
    right_eye_figs = [f for f in fig_files if f.is_file() and "RightEye" in f.name]
    
    print(f"Found {len(left_eye_figs)} figure files with LeftEye")
    print(f"Found {len(right_eye_figs)} figure files with RightEye")
    
    if len(left_eye_figs) > 0 or len(right_eye_figs) > 0:
        # Backup original figures directory
        if figures_backup_dir.exists():
            print(f"Warning: {figures_backup_dir} already exists. Skipping backup.")
        else:
            print("Creating backup of figures...")
            figures_dir.rename(figures_backup_dir)
            print(f"✓ Backed up to: {figures_backup_dir}")
        
        # Create new figures directory and swap files
        print("Swapping eye labels in figures...")
        figures_dir.mkdir(parents=True, exist_ok=True)
        
        for old_file in figures_backup_dir.glob("*"):
            if not old_file.is_file():
                continue
            
            filename = old_file.name
            
            if "LeftEye" in filename:
                # LeftEye -> RightEye
                newname = filename.replace("LeftEye", "RightEye")
                print(f"  {filename} -> {newname}")
                
            elif "RightEye" in filename:
                # RightEye -> LeftEye
                newname = filename.replace("RightEye", "LeftEye")
                print(f"  {filename} -> {newname}")
                
            else:
                # No eye label, copy as-is
                newname = filename
            
            # Use rsync for faster copying with progress
            src = str(old_file)
            dst = str(figures_dir / newname)
            os.system(f"rsync -av --progress {src} {dst}")
        
        print(f"✓ Figures processed successfully")
    else:
        print("No LeftEye/RightEye files found in figures directory")
else:
    print()
    print(f"Note: Figures directory not found at {figures_dir}")

print()
print("✓ Done! Eye labels swapped successfully.")
print()
print("Summary:")
print(f"  Original func files: {backup_dir}")
print(f"  New func files with swapped labels: {func_dir}")
if figures_dir.exists():
    print(f"  Original figures: {figures_backup_dir}")
    print(f"  New figures with swapped labels: {figures_dir}")

# Time
end_time = datetime.datetime.now()
print("\nStart time:\t{}\nEnd time:\t{}\nDuration:\t{}".format(
    start_time,
    end_time,
    end_time - start_time))

# Define permission cmd
print(f'\nChanging files permissions in {main_dir}/{project_dir}')
os.system(f"chmod -Rf 771 {main_dir}/{project_dir}")
os.system(f"chgrp -Rf {group} {main_dir}/{project_dir}")