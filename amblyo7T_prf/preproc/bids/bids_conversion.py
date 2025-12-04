"""
-----------------------------------------------------------------------------------------
bids_conversion.py
-----------------------------------------------------------------------------------------
Goal of the script:
Convert amblyopia 7T pRF data to BIDS format using automatic task mapping from Readme.txt
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject BIDS code (e.g. sub-03)
sys.argv[4]: subject source folder (e.g. bm)
sys.argv[5]: group (e.g. 327)
sys.argv[6]: subject group (A for amblyopes, C for controls)
-----------------------------------------------------------------------------------------
Output(s):
BIDS-formatted dataset with automatic task name mapping
Verification files in sourcedata folder
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/bids/
2. run python command
python bids_conversion.py [main directory] [project name] [subject BIDS] [subject code] [group] [subject group]
-----------------------------------------------------------------------------------------
Example:
python bids_conversion.py /scratch/mszinte/data amblyo7T_prf sub-01 bm 327 A
python bids_conversion.py /scratch/mszinte/data amblyo7T_prf sub-02 cg 327 A
python bids_conversion.py /scratch/mszinte/data amblyo7T_prf sub-03 hb 327 C
python bids_conversion.py /scratch/mszinte/data amblyo7T_prf sub-04 ii 327 C
python bids_conversion.py /scratch/mszinte/data amblyo7T_prf sub-05 kg 327 A
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# import modules
import sys
import os
import subprocess
import ipdb
deb = ipdb.set_trace
opj = os.path.join

# import utility functions
sys.path.append("{}/../../utils".format(os.getcwd()))
from bids_utils import (
    parse_readme,
    create_verification_file,
    run_dcm2niix,
    generate_psir_images,
    copy_anatomical,
    copy_functional_runs,
    copy_fieldmaps,
    update_functional_metadata,
    fix_nifti_tr,
    update_fieldmap_metadata,
    update_dataset_description
)

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
subject_code = sys.argv[4]
group = sys.argv[5]
subject_group = sys.argv[6]  # A or C

# define paths
base_path = opj(main_dir, project_dir)
sourcedata_dir = opj(base_path, 'sourcedata', subject_group, subject_code)
anatomy_dir = opj(base_path, 'sourcedata', f'{subject_group}-anatomy')
readme_path = opj(sourcedata_dir, 'Readme.txt')
dcm2niix_dir = opj(sourcedata_dir, 'dcm2niix')
parrec_dir = opj(sourcedata_dir, 'PARREC')
subject_dir = opj(base_path, subject)
session_dir = opj(subject_dir, 'ses-01')

# create output directories
os.makedirs(opj(session_dir, 'anat'), exist_ok=True)
os.makedirs(opj(session_dir, 'func'), exist_ok=True)
os.makedirs(opj(session_dir, 'fmap'), exist_ok=True)

print("=" * 80)
print("BIDS CONVERSION - AMBLYOPIA 7T pRF DATA")
print("=" * 80)
print(f"Subject: {subject}")
print(f"Source folder: {subject_code}")
print(f"Base path: {base_path}")
print()

# Step 1: Convert PARREC to NIfTI
print("Step 1: Converting PARREC to NIfTI...")
converted = run_dcm2niix(parrec_dir, dcm2niix_dir)
if converted:
    print("✓ Conversion complete\n")
else:
    print("✓ Skipped (already exists)\n")

# Step 2: Generate PSIR reconstructed images
print("Step 2: Generating PSIR reconstructed images (T1w)...")
try:
    psir_outputs = generate_psir_images(dcm2niix_dir)
    print("✓ PSIR images generated\n")
except FileNotFoundError as e:
    print(f"  WARNING: {e}")
    print("  Skipping PSIR generation\n")

# Step 3: Parse Readme.txt to get task mapping
print("Step 3: Parsing Readme.txt for task mapping...")
scans = parse_readme(readme_path)
print(f"✓ Found {len(scans)} scans\n")

# Step 4: Create verification file
print("Step 4: Creating verification file...")
verification_file = opj(sourcedata_dir, 'bids_mapping_verification.txt')
create_verification_file(scans, subject, readme_path, verification_file)
print(f"✓ Verification file created: {verification_file}\n")

# Display mapping
print("=" * 80)
print("BIDS CONVERSION MAPPING")
print("=" * 80)
with open(verification_file, 'r') as f:
    print(f.read())

# Ask for confirmation
response = input("\nDoes this mapping look correct? (y/n): ")
if response.lower() != 'y':
    print("Conversion aborted. Please check the Readme.txt file.")
    sys.exit(1)

# Step 5: Copy anatomical data
print("\nStep 5: Copying anatomical data...")
copy_anatomical(dcm2niix_dir, session_dir, subject)
print("✓ Anatomical data copied\n")

# Step 6: Extract and copy FreeSurfer mri folder
print("Step 6: Extracting and copying FreeSurfer mri folder...")
freesurfer_tgz = opj(anatomy_dir, f'{subject_code}.tgz')
extracted_folder = opj(anatomy_dir, subject_code)
freesurfer_dir = opj(base_path, 'derivatives', 'fmriprep', 'freesurfer', f'{subject}_ses-01')

if os.path.exists(freesurfer_tgz):
    # Unzip the tgz file if extracted folder doesn't exist
    if not os.path.exists(extracted_folder):
        print(f"  Extracting {freesurfer_tgz}...")
        subprocess.run(['tar', '-xzf', freesurfer_tgz, '-C', anatomy_dir], check=True)
        print("  ✓ Extraction complete")
    else:
        print("  ✓ FreeSurfer folder already extracted")
    
    # Copy contents of extracted folder to FreeSurfer derivatives
    if os.path.exists(extracted_folder):
        os.makedirs(freesurfer_dir, exist_ok=True)
        # Copy contents with cp -r source/. dest/
        subprocess.run(['cp', '-r', f'{extracted_folder}/.', freesurfer_dir], check=True)
        print(f"  ✓ FreeSurfer files copied to {freesurfer_dir}\n")
    else:
        print(f"  WARNING: Extracted folder not found at {extracted_folder}\n")
else:
    print(f"  WARNING: FreeSurfer archive not found at {freesurfer_tgz}\n")

# Step 7: Copy functional runs
print("Step 7: Copying functional runs...")
copy_commands_file = opj(sourcedata_dir, 'copy_commands.sh')
copy_functional_runs(scans, dcm2niix_dir, session_dir, subject, copy_commands_file)
print("✓ Functional runs copied\n")

# Step 8: Copy fieldmaps
print("Step 8: Copying fieldmaps...")
copy_fieldmaps(dcm2niix_dir, session_dir, subject)
print("✓ Fieldmaps copied\n")

# Step 9: Update functional metadata
print("Step 9: Updating functional JSON metadata...")
update_functional_metadata(session_dir, subject)
print("✓ Functional metadata updated\n")

# Step 10: Fix TR in NIfTI headers
print("Step 10: Fixing TR in NIfTI headers...")
fix_nifti_tr(session_dir, subject)
print("✓ TR fixed in NIfTI headers\n")

# Step 11: Update fieldmap metadata
print("Step 11: Updating fieldmap metadata...")
update_fieldmap_metadata(session_dir, subject)
print("✓ Fieldmap metadata updated\n")

# Step 12: Update dataset description
print("Step 12: Updating dataset_description.json...")
update_dataset_description(base_path)
print("✓ Dataset description updated\n")

# Step 13: Set permissions
print("Step 13: Setting permissions...")
os.system(f"chmod -Rf 771 {base_path}")
os.system(f"chgrp -Rf {group} {base_path}")
print("✓ Permissions set\n")

print("=" * 80)
print("BIDS CONVERSION COMPLETE")
print("=" * 80)
print(f"Mapping verification file: {verification_file}")
print("Please verify the output with BIDS validator")
