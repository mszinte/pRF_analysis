#!/bin/bash

# Subject ID variable
SUBJECT="sub-01"

# Copy manually from local to sourcedata

# Permissions
chmod -Rf 771 /scratch/mszinte/data/amblyo7T_prf/
chgrp -Rf 327 /scratch/mszinte/data/amblyo7T_prf/

# Convert PARREC to nifti
cd /scratch/mszinte/data/amblyo7T_prf/
mkdir sourcedata/BM/dcm2niix/
dcm2niix -o sourcedata/BM/dcm2niix/ -z y sourcedata/BM/PARREC/


# create manually readme
# create manually participants.json and .tsv
# create manually dataset_description.json

# ${SUBJECT}
# create subject folder
# create session folder
mkdir ${SUBJECT}/
mkdir ${SUBJECT}/ses-01/
mkdir ${SUBJECT}/ses-01/anat/
mkdir ${SUBJECT}/ses-01/fmap/
mkdir ${SUBJECT}/ses-01/func/


# Define source and target directories
cd /scratch/mszinte/data/amblyo7T_prf/
SOURCE_DIR="sourcedata/BM/dcm2niix"
SUB_DIR="${SUBJECT}/ses-01"

# Copy anatomical (T1w PSIR second inversion)
cp ${SOURCE_DIR}/PARREC_WIPPSIR1mm3SENSE_20141016132749_15_t2387.nii.gz ${SUB_DIR}/anat/${SUBJECT}_ses-01_T1w.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPPSIR1mm3SENSE_20141016132749_15_t2387.json ${SUB_DIR}/anat/${SUBJECT}_ses-01_T1w.json 

# Copy functional runs
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_5.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-WedgesLeftEye_run-01_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_6.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-01_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_7.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-RingsLeftEye_run-01_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_8.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-02_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_10.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsRightEye_run-01_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_11.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-WedgesRightEye_run-01_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_12.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-RingsRightEye_run-01_bold.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_13.nii.gz ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsRightEye_run-02_bold.nii.gz

# Copy functional JSON sidecars
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_5.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-WedgesLeftEye_run-01_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_6.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-01_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_7.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-RingsLeftEye_run-01_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_8.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-02_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_10.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsRightEye_run-01_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_11.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-WedgesRightEye_run-01_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_12.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-RingsRightEye_run-01_bold.json
cp ${SOURCE_DIR}/PARREC_WIPfMRI_RL_33_166SENSE_20141016132749_13.json ${SUB_DIR}/func/${SUBJECT}_ses-01_task-BarsRightEye_run-02_bold.json

# Copy fieldmap files
cp ${SOURCE_DIR}/PARREC_WIPt2s_phase_0p548_slicesSENSE_20141016132749_14_ph.nii.gz ${SUB_DIR}/fmap/${SUBJECT}_ses-01_fieldmap.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPt2s_phase_0p548_slicesSENSE_20141016132749_14_ph.json ${SUB_DIR}/fmap/${SUBJECT}_ses-01_fieldmap.json
cp ${SOURCE_DIR}/PARREC_WIPt2s_phase_0p548_slicesSENSE_20141016132749_14.nii.gz ${SUB_DIR}/fmap/${SUBJECT}_ses-01_magnitude.nii.gz


## CORRECTION AFTER CHECKING BIDS VALIDATOR

# Add TaskName and RepetitionTime to each JSON file
python3 << EOF
import json

tasks = [
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-WedgesLeftEye_run-01_bold.json", "WedgesLeftEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-01_bold.json", "BarsLeftEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-RingsLeftEye_run-01_bold.json", "RingsLeftEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-02_bold.json", "BarsLeftEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsRightEye_run-01_bold.json", "BarsRightEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-WedgesRightEye_run-01_bold.json", "WedgesRightEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-RingsRightEye_run-01_bold.json", "RingsRightEye"),
    ("${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsRightEye_run-02_bold.json", "BarsRightEye")
]

for json_file, task_name in tasks:
    with open(json_file, 'r') as f:
        data = json.load(f)
    data['TaskName'] = task_name
    data['RepetitionTime'] = 2.0
    data['SkullStripped'] = True
    data['PhaseEncodingDirection'] = 'i-'
    data['Manufacturer'] = 'Philips'
    data['ManufacturersModelName'] = 'Intera Achieva'
    data['MagneticFieldStrength'] = 7.0
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)
EOF


# Fix TR in NIfTI headers using Python
python3 << EOF
import nibabel as nib
import numpy as np

func_files = [
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-WedgesLeftEye_run-01_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-01_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-RingsLeftEye_run-01_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsLeftEye_run-02_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsRightEye_run-01_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-WedgesRightEye_run-01_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-RingsRightEye_run-01_bold.nii.gz",
    "${SUBJECT}/ses-01/func/${SUBJECT}_ses-01_task-BarsRightEye_run-02_bold.nii.gz"
]

for func_file in func_files:
    img = nib.load(func_file)
    hdr = img.header.copy()
    zooms = list(hdr.get_zooms())
    zooms[3] = 2.0
    hdr.set_zooms(zooms)
    new_img = nib.Nifti1Image(img.get_fdata(), img.affine, hdr)
    nib.save(new_img, func_file)
EOF

# Add Units to fieldmap JSON
python3 << EOF
import json
import os
import glob

# Get all functional bold.nii.gz files
func_dir = '${SUBJECT}/ses-01/func'
func_files = sorted(glob.glob(os.path.join(func_dir, '*_bold.nii.gz')))

# Create IntendedFor list with relative paths from subject directory
intended_for = []
for func_file in func_files:
    # Get relative path from subject directory (remove subject prefix)
    rel_path = func_file.replace('${SUBJECT}/', '')
    intended_for.append(rel_path)

with open('${SUBJECT}/ses-01/fmap/${SUBJECT}_ses-01_fieldmap.json', 'r') as f:
    data = json.load(f)
data['Units'] = 'Hz'
data['SkullStripped'] = False
data['Manufacturer'] = 'Philips'
data['ManufacturersModelName'] = 'Intera Achieva'
data['MagneticFieldStrength'] = 7.0
data['IntendedFor'] = intended_for
with open('${SUBJECT}/ses-01/fmap/${SUBJECT}_ses-01_fieldmap.json', 'w') as f:
    json.dump(data, f, indent=2)

# Update T1w JSON sidecar
with open('${SUBJECT}/ses-01/anat/${SUBJECT}_ses-01_T1w.json', 'r') as f:
    t1w_data = json.load(f)
t1w_data['SkullStripped'] = False
t1w_data['Manufacturer'] = 'Philips'
t1w_data['ManufacturersModelName'] = 'Intera Achieva'
t1w_data['MagneticFieldStrength'] = 7.0
with open('${SUBJECT}/ses-01/anat/${SUBJECT}_ses-01_T1w.json', 'w') as f:
    json.dump(t1w_data, f, indent=2)
EOF


# Add License to dataset_description.json
python3 << 'EOF'
import json

try:
    with open('dataset_description.json', 'r') as f:
        data = json.load(f)
except FileNotFoundError:
    data = {}

data['License'] = 'CC0'
data['GeneratedBy'] = [
    {
        "Name": "Schluppeck, D., Arnoldussen, D., Hussain, Z., Besle, J., Francis, S. T., & McGraw, P. V.",
        "Description": "Strabismus and amblyopia disrupt spatial perception but not the fidelity of cortical maps in human primary visual cortex"
    }
]
data['SourceDatasets'] = [
    {
        "DOI": "10.1016/j.visres.2025.108677",
        "URL": "https://doi.org/10.1016/j.visres.2025.108677"
    }
]

with open('dataset_description.json', 'w') as f:
    json.dump(data, f, indent=2)
EOF


# Permissions
chmod -Rf 771 /scratch/mszinte/data/amblyo7T_prf/
chgrp -Rf 327 /scratch/mszinte/data/amblyo7T_prf/
