#!/bin/bash

# Subject ID variable
SUBJECT="sub-03"

# Path to Readme.txt (should be in sourcedata folder before PARREC/)
README_PATH="/scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/Readme.txt"

# Copy manually from local to sourcedata

# Permissions
chmod -Rf 771 /scratch/mszinte/data/amblyo7T_prf/
chgrp -Rf 327 /scratch/mszinte/data/amblyo7T_prf/

# Convert PARREC to nifti
cd /scratch/mszinte/data/amblyo7T_prf/
mkdir -p sourcedata/BM/dcm2niix/
dcm2niix -o sourcedata/BM/dcm2niix/ -z y sourcedata/BM/PARREC/

# Create T1images
# RUN ../anatomical/psir_mougin_adapted.ipynb

# create manually readme
# create manually participants.json and .tsv
# create manually dataset_description.json

# Parse Readme.txt and create mapping file
python3 << EOF
import re
import json

readme_path = "${README_PATH}"
subject = "${SUBJECT}"

# Read and parse Readme.txt
with open(readme_path, 'r') as f:
    lines = f.readlines()

# Find the table start
table_start = False
scans = []

for line in lines:
    if line.strip().startswith('scan\tdescription'):
        table_start = True
        continue
    
    if table_start and line.strip():
        parts = line.strip().split('\t')
        if len(parts) >= 3 and parts[0].isdigit():
            scan_num = int(parts[0])
            description = parts[1].strip().upper()
            filename = parts[2].strip()
            
            # Extract file number from filename (e.g., ViasAmb77_5_1 -> 5)
            file_match = re.search(r'_(\d+)_\d+_modulus', filename)
            if file_match:
                file_num = file_match.group(1)
            else:
                file_num = None
            
            # Parse description to get task and eye
            desc_parts = description.replace('_', ' ').split()
            if len(desc_parts) >= 2:
                task_type = desc_parts[0]  # RINGS, BARS, WEDGES, etc.
                eye = desc_parts[1]  # LEFT or RIGHT
                
                # Create BIDS task name: TaskTypeEye format (e.g., WedgesLeftEye)
                task_name = f"{task_type.capitalize()}{eye.capitalize()}Eye"
                
                scans.append({
                    'scan_num': scan_num,
                    'file_num': file_num,
                    'description': description,
                    'task_name': task_name,
                    'eye': eye,
                    'task_type': task_type
                })

# Count runs for each task
task_counts = {}
for scan in scans:
    task = scan['task_name']
    if task not in task_counts:
        task_counts[task] = 0
    task_counts[task] += 1
    scan['run_num'] = task_counts[task]

# Create human-readable mapping file in sourcedata/bm/
mapping_file = '/scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/bids_mapping_verification.txt'
with open(mapping_file, 'w') as f:
    f.write("BIDS CONVERSION MAPPING - VERIFICATION FILE\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Subject: {subject}\n")
    f.write(f"Readme file: {readme_path}\n\n")
    f.write(f"{'Scan':<6} {'File#':<8} {'Description':<20} {'BIDS Task Name':<30} {'Run':<5}\n")
    f.write("-" * 80 + "\n")
    
    for scan in scans:
        f.write(f"{scan['scan_num']:<6} "
                f"{scan['file_num']:<8} "
                f"{scan['description']:<20} "
                f"task-{scan['task_name']}_run-{scan['run_num']:02d}_bold"
                f"{'':<5}\n")
    
    f.write("\n" + "=" * 80 + "\n")

print(f"Mapping file created: {mapping_file}")

# Save scans data for later use
import pickle
with open('/tmp/bids_scans_mapping.pkl', 'wb') as f:
    pickle.dump(scans, f)

EOF

# Display mapping for verification
echo ""
echo "=========================================="
echo "BIDS CONVERSION MAPPING"
echo "=========================================="
cat /scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/bids_mapping_verification.txt
echo ""
read -p "Does this mapping look correct? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Conversion aborted. Please check the Readme.txt file."
    exit 1
fi

# ${SUBJECT}
# create subject folder
cd /scratch/mszinte/data/amblyo7T_prf/
mkdir -p ${SUBJECT}/
mkdir -p ${SUBJECT}/ses-01/
mkdir -p ${SUBJECT}/ses-01/anat/
mkdir -p ${SUBJECT}/ses-01/fmap/
mkdir -p ${SUBJECT}/ses-01/func/

# Define source and target directories
cd /scratch/mszinte/data/amblyo7T_prf/
SOURCE_DIR="sourcedata/BM/dcm2niix"
SUB_DIR="${SUBJECT}/ses-01"

# Copy anatomical (PSIRbeta_Mougin out of notebook - see above)
cp ${SOURCE_DIR}/PARREC_WIPPSIR1mm3SENSE_20141016132749_15_PSIRbeta_Mougin.nii.gz ${SUB_DIR}/anat/${SUBJECT}_ses-01_T1w.nii.gz
cp ${SOURCE_DIR}/PARREC_WIPPSIR1mm3SENSE_20141016132749_15_t2387.json ${SUB_DIR}/anat/${SUBJECT}_ses-01_T1w.json 

# Copy functional runs and JSON sidecars using mapping
python3 << EOF
import pickle
import os
import subprocess

with open('/tmp/bids_scans_mapping.pkl', 'rb') as f:
    scans = pickle.load(f)

source_dir = "${SOURCE_DIR}"
sub_dir = "${SUB_DIR}"
subject = "${SUBJECT}"

# Find the base pattern for dcm2niix output files
# This needs to be adapted based on actual dcm2niix output naming
# Assuming pattern: PARREC_WIPfMRI_RL_33_166SENSE_YYYYMMDDHHMMSS_{file_num}.nii.gz

# First, list all files in source directory to find the pattern
dcm2niix_files = {}
for filename in os.listdir(source_dir):
    if filename.endswith('.nii.gz') and 'fMRI' in filename:
        # Extract the file number from the end of filename
        parts = filename.replace('.nii.gz', '').split('_')
        if parts[-1].isdigit():
            file_num = parts[-1]
            dcm2niix_files[file_num] = filename.replace('.nii.gz', '')

print("Found dcm2niix files:")
for num, base in sorted(dcm2niix_files.items()):
    print(f"  File #{num}: {base}")

print("\nCopying functional files...")
copy_commands = []

for scan in scans:
    file_num = scan['file_num']
    task_name = scan['task_name']
    run_num = scan['run_num']
    
    # Try both with and without leading zero
    file_num_int = str(int(file_num))  # Remove leading zeros
    
    source_base = None
    if file_num in dcm2niix_files:
        source_base = dcm2niix_files[file_num]
    elif file_num_int in dcm2niix_files:
        source_base = dcm2niix_files[file_num_int]
    
    if source_base:
        target_base = f"{subject}_ses-01_task-{task_name}_run-{run_num:02d}_bold"
        
        # NIfTI file
        src_nii = f"{source_dir}/{source_base}.nii.gz"
        dst_nii = f"{sub_dir}/func/{target_base}.nii.gz"
        copy_commands.append(f"cp {src_nii} {dst_nii}")
        
        # JSON file
        src_json = f"{source_dir}/{source_base}.json"
        dst_json = f"{sub_dir}/func/{target_base}.json"
        copy_commands.append(f"cp {src_json} {dst_json}")
        
        print(f"  Scan {scan['scan_num']}: {source_base} -> {target_base}")
    else:
        print(f"  WARNING: File #{file_num} not found for scan {scan['scan_num']}")

# Write copy commands to shell script
with open('/scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/copy_commands.sh', 'w') as f:
    f.write('#!/bin/bash\n\n')
    for cmd in copy_commands:
        f.write(cmd + '\n')

print("\nCopy commands written to /scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/copy_commands.sh")

EOF

# Execute copy commands
chmod +x /scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/copy_commands.sh
/scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/copy_commands.sh

# Copy fieldmap files
cp ${SOURCE_DIR}/PARREC_WIPt2s_phase_0p548_slicesSENSE_*_14_ph.nii.gz ${SUB_DIR}/fmap/${SUBJECT}_ses-01_fieldmap.nii.gz 2>/dev/null || echo "Warning: Fieldmap phase not found"
cp ${SOURCE_DIR}/PARREC_WIPt2s_phase_0p548_slicesSENSE_*_14_ph.json ${SUB_DIR}/fmap/${SUBJECT}_ses-01_fieldmap.json 2>/dev/null || echo "Warning: Fieldmap phase JSON not found"
cp ${SOURCE_DIR}/PARREC_WIPt2s_phase_0p548_slicesSENSE_*_14.nii.gz ${SUB_DIR}/fmap/${SUBJECT}_ses-01_magnitude.nii.gz 2>/dev/null || echo "Warning: Fieldmap magnitude not found"


## CORRECTION AFTER CHECKING BIDS VALIDATOR

# Add TaskName and other metadata to each JSON file
python3 << EOF
import json
import glob
import os

# Get all functional bold JSON files
func_jsons = glob.glob('${SUBJECT}/ses-01/func/*_bold.json')

for json_file in func_jsons:
    # Extract task name from filename
    basename = os.path.basename(json_file)
    task_match = basename.split('task-')[1].split('_')[0]
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    data['TaskName'] = task_match
    data['RepetitionTime'] = 2.0
    data['SkullStripped'] = True
    data['PhaseEncodingDirection'] = 'i-'
    data['Manufacturer'] = 'Philips'
    data['ManufacturersModelName'] = 'Intera Achieva'
    data['MagneticFieldStrength'] = 7.0
    
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)

print(f"Updated {len(func_jsons)} functional JSON files")
EOF


# Fix TR in NIfTI headers using Python
python3 << EOF
import nibabel as nib
import glob

func_files = glob.glob('${SUBJECT}/ses-01/func/*_bold.nii.gz')

for func_file in func_files:
    img = nib.load(func_file)
    hdr = img.header.copy()
    zooms = list(hdr.get_zooms())
    zooms[3] = 2.0
    hdr.set_zooms(zooms)
    new_img = nib.Nifti1Image(img.get_fdata(), img.affine, hdr)
    nib.save(new_img, func_file)

print(f"Fixed TR in {len(func_files)} NIfTI files")
EOF

# Add Units to fieldmap JSON and create IntendedFor
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

# Update fieldmap JSON
fieldmap_json = '${SUBJECT}/ses-01/fmap/${SUBJECT}_ses-01_fieldmap.json'
if os.path.exists(fieldmap_json):
    with open(fieldmap_json, 'r') as f:
        data = json.load(f)
    data['Units'] = 'Hz'
    data['SkullStripped'] = False
    data['Manufacturer'] = 'Philips'
    data['ManufacturersModelName'] = 'Intera Achieva'
    data['MagneticFieldStrength'] = 7.0
    data['IntendedFor'] = intended_for
    with open(fieldmap_json, 'w') as f:
        json.dump(data, f, indent=2)
    print("Updated fieldmap JSON")

# Update T1w JSON sidecar
t1w_json = '${SUBJECT}/ses-01/anat/${SUBJECT}_ses-01_T1w.json'
if os.path.exists(t1w_json):
    with open(t1w_json, 'r') as f:
        t1w_data = json.load(f)
    t1w_data['SkullStripped'] = False
    t1w_data['Manufacturer'] = 'Philips'
    t1w_data['ManufacturersModelName'] = 'Intera Achieva'
    t1w_data['MagneticFieldStrength'] = 7.0
    with open(t1w_json, 'w') as f:
        json.dump(t1w_data, f, indent=2)
    print("Updated T1w JSON")
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

echo ""
echo "=========================================="
echo "BIDS CONVERSION COMPLETE"
echo "=========================================="
echo "Mapping verification file: /scratch/mszinte/data/amblyo7T_prf/sourcedata/bm/bids_mapping_verification.txt"
echo "Please verify the output with BIDS validator"
