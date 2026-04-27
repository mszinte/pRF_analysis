"""
-----------------------------------------------------------------------------------------
bids_conversion_special.py
-----------------------------------------------------------------------------------------
Goal of the script:
Convert amblyopia 7T pRF data to BIDS format for special subjects (sub-05, sub-15)
that differ from the standard pipeline in the following ways:
  - Raw NIfTI files already exist (no PARREC conversion needed)
  - No JSON sidecar files -> copied from sub-01 BIDS output and adapted
  - NIfTI files are not gzipped -> gzipped in-place
  - No fieldmaps available
  - No anatomical data -> T1w derived from FreeSurfer T1.mgz
  - sub-15 has only 6 functional runs instead of 8
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject BIDS code (e.g. sub-05)
sys.argv[4]: subject source folder (e.g. sm)
sys.argv[5]: group (e.g. 327)
sys.argv[6]: subject group (A for amblyopes, C for controls)
-----------------------------------------------------------------------------------------
Output(s):
BIDS-formatted dataset
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/bids/
2. run python command
python bids_conversion_special.py [main directory] [project name] [subject BIDS] [subject code] [group] [subject group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/bids/
python bids_conversion_special.py /scratch/mszinte/data amblyo7T_prf sub-05 xx 327 A
python bids_conversion_special.py /scratch/mszinte/data amblyo7T_prf sub-15 xx 327 C
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

import sys
import os
import json
import subprocess
import ipdb
deb = ipdb.set_trace
opj = os.path.join

sys.path.append("{}/../../utils".format(os.getcwd()))
from bids_utils import (
    parse_readme,
    create_verification_file,
    convert_freesurfer_t1_to_nifti,
    update_dataset_description,
    gzip_nifti_files,
    copy_and_adapt_anat_json,
    copy_functional_runs_special,
    fix_nifti_tr,
)

# inputs
main_dir      = sys.argv[1]
project_dir   = sys.argv[2]
subject       = sys.argv[3]   # e.g. sub-05
subject_code  = sys.argv[4]   # e.g. sm
group         = sys.argv[5]   # e.g. 327
subject_group = sys.argv[6]   # A or C

# define paths
base_path      = opj(main_dir, project_dir)
sourcedata_dir = opj(base_path, 'sourcedata', subject_group, subject_code)
anatomy_dir    = opj(base_path, 'sourcedata', f'{subject_group}-anatomy')
readme_path    = opj(sourcedata_dir, 'Readme.txt')
raw_dir        = opj(sourcedata_dir, 'Raw')
subject_dir    = opj(base_path, subject)
session_dir    = opj(subject_dir, 'ses-01')
freesurfer_dir = opj(base_path, 'derivatives', 'fmriprep', 'freesurfer', f'{subject}_ses-01')

# sub-01 template directories for JSON sidecars
sub01_func_dir = opj(base_path, 'sub-01', 'ses-01', 'func')
sub01_anat_dir = opj(base_path, 'sub-01', 'ses-01', 'anat')

# create BIDS output directories (no fmap for these subjects)
os.makedirs(opj(session_dir, 'anat'), exist_ok=True)
os.makedirs(opj(session_dir, 'func'), exist_ok=True)

print("=" * 80)
print("BIDS CONVERSION (SPECIAL) - AMBLYOPIA 7T pRF DATA")
print("=" * 80)
print(f"Subject      : {subject}")
print(f"Source folder: {subject_code}")
print(f"Raw NIfTI dir: {raw_dir}")
print(f"Base path    : {base_path}")
print()

# --------------------------------------------------------------------------
# Step 1: Gzip raw NIfTI files
# --------------------------------------------------------------------------
print("Step 1: Gzipping raw NIfTI files...")
gz_files = gzip_nifti_files(raw_dir)
print(f"✓ {len(gz_files)} .nii.gz files ready\n")

# --------------------------------------------------------------------------
# Step 2: Parse Readme.txt
# --------------------------------------------------------------------------
print("Step 2: Parsing Readme.txt for task mapping...")
scans = parse_readme(readme_path)
print(f"✓ Found {len(scans)} scans\n")

# --------------------------------------------------------------------------
# Step 3: Create and display verification file
# --------------------------------------------------------------------------
print("Step 3: Creating verification file...")
verification_file = opj(sourcedata_dir, 'bids_mapping_verification.txt')
create_verification_file(scans, subject, readme_path, verification_file)
print(f"✓ Verification file created: {verification_file}\n")

print("=" * 80)
print("BIDS CONVERSION MAPPING")
print("=" * 80)
with open(verification_file, 'r') as f:
    print(f.read())

response = input("\nDoes this mapping look correct? (y/n): ")
if response.lower() != 'y':
    print("Conversion aborted. Please check the Readme.txt file.")
    sys.exit(1)

# --------------------------------------------------------------------------
# Step 4: Extract FreeSurfer archive and convert T1.mgz -> T1w NIfTI
# --------------------------------------------------------------------------
print("\nStep 4: Extracting FreeSurfer archive...")
freesurfer_tgz   = opj(anatomy_dir, f'{subject_code}.tgz')
extracted_folder = opj(anatomy_dir, subject_code)

if os.path.exists(freesurfer_tgz):
    if not os.path.exists(extracted_folder):
        print(f"  Extracting {freesurfer_tgz}...")
        subprocess.run(['tar', '-xzf', freesurfer_tgz, '-C', anatomy_dir], check=True)
        print("  ✓ Extraction complete")
    else:
        print("  ✓ FreeSurfer folder already extracted")

    os.makedirs(freesurfer_dir, exist_ok=True)
    subprocess.run(['cp', '-r', f'{extracted_folder}/.', freesurfer_dir], check=True)
    print(f"  ✓ FreeSurfer files copied to {freesurfer_dir}\n")
else:
    print(f"  WARNING: FreeSurfer archive not found at {freesurfer_tgz}\n")

print("Step 4b: Converting FreeSurfer T1.mgz to NIfTI T1w...")
t1w_path = opj(session_dir, 'anat', f'{subject}_ses-01_T1w.nii.gz')
if not os.path.exists(t1w_path):
    convert_freesurfer_t1_to_nifti(freesurfer_dir, session_dir, subject)
else:
    print("  ✓ T1w already exists, skipping\n")

# Copy anat JSON from sub-01 template and flag the source
anat_json_path = opj(session_dir, 'anat', f'{subject}_ses-01_T1w.json')
if not os.path.exists(anat_json_path):
    print("  Copying T1w JSON template from sub-01...")
    try:
        copy_and_adapt_anat_json(sub01_anat_dir, anat_json_path)
        with open(anat_json_path, 'r') as f:
            anat_meta = json.load(f)
        anat_meta['Source'] = 'FreeSurfer T1.mgz'
        with open(anat_json_path, 'w') as f:
            json.dump(anat_meta, f, indent=2)
        print("  ✓ T1w JSON created\n")
    except FileNotFoundError as e:
        print(f"  WARNING: {e}\n")
else:
    print("  ✓ T1w JSON already exists\n")

# --------------------------------------------------------------------------
# Step 5: Copy functional runs + generate JSON sidecars
# --------------------------------------------------------------------------
print("Step 5: Copying functional runs and generating JSON sidecars...")
copy_commands_file = opj(sourcedata_dir, 'copy_commands.sh')
copy_functional_runs_special(scans, raw_dir, session_dir, subject,
                              sub01_func_dir, copy_commands_file)
print("✓ Functional runs copied\n")

# --------------------------------------------------------------------------
# Step 6: Fix TR in NIfTI headers
# --------------------------------------------------------------------------
print("Step 6: Fixing TR in NIfTI headers...")
fix_nifti_tr(session_dir, subject)
print("✓ TR fixed\n")

# --------------------------------------------------------------------------
# Step 7: Update dataset_description.json
# --------------------------------------------------------------------------
print("Step 7: Updating dataset_description.json...")
update_dataset_description(base_path)
print("✓ Dataset description updated\n")

# --------------------------------------------------------------------------
# Step 8: Set permissions
# --------------------------------------------------------------------------
print("Step 8: Setting permissions...")
os.system(f"chmod -Rf 771 {base_path}")
os.system(f"chgrp -Rf {group} {base_path}")
print("✓ Permissions set\n")

print("=" * 80)
print("BIDS CONVERSION (SPECIAL) COMPLETE")
print("=" * 80)
print(f"Subject              : {subject}")
print(f"Functional runs      : {len(scans)}")
print(f"Fieldmaps            : none (not available for this subject)")
print(f"Mapping verification : {verification_file}")
print("Please verify the output with BIDS validator")