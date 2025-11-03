#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:48:39 2023
Modified: Mon Oct 27 2025
-----------------------------------------------------------------------------------------
organise_presurfer.py
-----------------------------------------------------------------------------------------
Goal of the script:
Copy and rename presurfer T1w output to BIDS format
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: mesocentre project data directory
sys.argv[2]: individual subject name (e.g., sub-t044)
-----------------------------------------------------------------------------------------
To run:
On mesocentre
>> python organise_presurfer.py [meso_proj_dir] [subject]
python organise_presurfer.py /scratch/mszinte/data/centbids sub-t044
-----------------------------------------------------------------------------------------
written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""
import os
import sys
import shutil

# Define inputs
project_dir = sys.argv[1]
subject = sys.argv[2]
sub_num = subject.split('-')[-1]

print(f"Processing subject: {sub_num}")

# Check project directory exists
if not os.path.isdir(project_dir):
    print(f"✗ Project directory doesn't exist: {project_dir}")
    sys.exit(1)

print("✓ Project directory found")

# Define directories
presurfer_dir = os.path.join(project_dir, "derivatives/presurfer/presurfer")
subject_dir = os.path.join(presurfer_dir, f"sub-{sub_num}")

# Check subject directory exists
if not os.path.isdir(subject_dir):
    print(f"✗ Subject directory not found: {subject_dir}")
    sys.exit(1)

# Get session
sessions = os.listdir(subject_dir)
if not sessions:
    print(f"✗ No sessions found for subject {sub_num}")
    sys.exit(1)

session = sessions[0]
print(f"Session: {session}")

# Define source and destination directories
source_anat_dir = os.path.join(subject_dir, session, "anat")
dest_anat_dir = os.path.join(project_dir, f"sub-{sub_num}", session, "anat")

# Check source directory exists
if not os.path.isdir(source_anat_dir):
    print(f"✗ Source anatomy directory not found: {source_anat_dir}")
    sys.exit(1)

# Create destination directory
os.makedirs(dest_anat_dir, exist_ok=True)
print(f"✓ Destination directory: {dest_anat_dir}")

# Find the unstripped MPRAGEised file
mprageised_files = [f for f in os.listdir(source_anat_dir) if "desc-unstripped_MPRAGEised" in f]

if not mprageised_files:
    print(f"✗ No presurfer output found in {source_anat_dir}")
    sys.exit(1)

# Copy and rename T1w file
source_file = os.path.join(source_anat_dir, mprageised_files[0])

# Create new filename: sub-{sub_num}_{session}_T1w.nii.gz
dest_filename = f"sub-{sub_num}_{session}_T1w.nii.gz"
dest_file = os.path.join(dest_anat_dir, dest_filename)

# Copy file
shutil.copy2(source_file, dest_file)
print(f"✓ Copied T1w: {dest_filename}")
print(f"✓ Source: {source_file}")
print(f"✓ Destination: {dest_file}")
print("✓ All done!")