import os
import re
import json
import glob
import sys
import subprocess
import nibabel as nib
import numpy as np

opj = os.path.join

def update_functional_metadata(session_dir, subject):
    """
    Update functional JSON files with required BIDS metadata.
    
    Parameters
    ----------
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    """
    func_jsons = glob.glob(opj(session_dir, 'func', '*_bold.json'))
    
    for json_file in func_jsons:
        # Extract task name from filename
        basename = os.path.basename(json_file)
        task_match = basename.split('task-')[1].split('_')[0]
        
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        data['TaskName'] = task_match
        data['RepetitionTime'] = 1.6
        data['SkullStripped'] = True
        data['PhaseEncodingDirection'] = 'j-'
        data['Manufacturer'] = 'Philips'
        data['MagneticFieldStrength'] = 7.0
        data['TotalReadoutTime'] = 0.071
        
        with open(json_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    print(f"  Updated {len(func_jsons)} functional JSON files")

def fix_nifti_tr(session_dir, subject):
    """
    Fix TR in NIfTI headers for functional files.
    
    Parameters
    ----------
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    """
    func_files = glob.glob(opj(session_dir, 'func', '*_bold.nii.gz'))
    
    for func_file in func_files:
        img = nib.load(func_file)
        hdr = img.header.copy()
        zooms = list(hdr.get_zooms())
        zooms[3] = 1.6  # Match the TR in JSON (1.6s)
        hdr.set_zooms(zooms)
        new_img = nib.Nifti1Image(img.get_fdata(), img.affine, hdr)
        nib.save(new_img, func_file)
    
    print(f"  Fixed TR in {len(func_files)} NIfTI files")

def update_epi_fieldmap_metadata(session_dir, subject):
    """
    Update EPI fieldmap JSON files with required BIDS metadata.
    These are the dir-AP and dir-PA epi files used for distortion correction.
    
    Parameters
    ----------
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    """
    # Get all functional files for IntendedFor field
    func_files = sorted(glob.glob(opj(session_dir, 'func', '*_bold.nii.gz')))
    
    # Create IntendedFor list with relative paths
    intended_for = []
    for func_file in func_files:
        rel_path = func_file.replace(f"{os.path.dirname(session_dir)}/", '')
        intended_for.append(rel_path)
    
    # Update all fieldmap epi JSON files (dir-AP and dir-PA)
    fmap_jsons = glob.glob(opj(session_dir, 'fmap', '*_epi.json'))
    
    for fmap_json in fmap_jsons:
        basename = os.path.basename(fmap_json)
        
        # Read existing data or create new dict
        if os.path.exists(fmap_json):
            with open(fmap_json, 'r') as f:
                data = json.load(f)
        else:
            data = {}
        
        # Set common metadata
        data['SkullStripped'] = True
        data['Manufacturer'] = 'Philips'
        data['MagneticFieldStrength'] = 7.0
        data['TotalReadoutTime'] = 0.071
        data['IntendedFor'] = intended_for
        
        # Set PhaseEncodingDirection based on filename
        if 'dir-AP' in basename:
            data['PhaseEncodingDirection'] = 'j'
        elif 'dir-PA' in basename:
            data['PhaseEncodingDirection'] = 'j-'
        
        with open(fmap_json, 'w') as f:
            json.dump(data, f, indent=2)
    
    print(f"  Updated {len(fmap_jsons)} EPI fieldmap JSON files")

def update_fieldmap_metadata(session_dir, subject):
    """
    Update legacy fieldmap JSON with required BIDS metadata and IntendedFor field.
    
    Parameters
    ----------
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    """
    # Get all functional files
    func_files = sorted(glob.glob(opj(session_dir, 'func', '*_bold.nii.gz')))
    
    # Create IntendedFor list with relative paths
    intended_for = []
    for func_file in func_files:
        rel_path = func_file.replace(f"{os.path.dirname(session_dir)}/", '')
        intended_for.append(rel_path)
    
    # Update old-style fieldmap JSON if it exists
    fieldmap_jsons = glob.glob(opj(session_dir, 'fmap', '*_fieldmap.json'))
    
    for fieldmap_json in fieldmap_jsons:
        with open(fieldmap_json, 'r') as f:
            data = json.load(f)
        
        data['Units'] = 'Hz'
        data['SkullStripped'] = True
        data['Manufacturer'] = 'Philips'
        data['MagneticFieldStrength'] = 7.0
        data['IntendedFor'] = intended_for
        
        with open(fieldmap_json, 'w') as f:
            json.dump(data, f, indent=2)
    
    if fieldmap_jsons:
        print(f"  Updated {len(fieldmap_jsons)} legacy fieldmap JSON files")
    # Step 5: Update anatomical metadata
    print("Step 5: Updating anatomical metadata...")
    
    # Update T1w JSON
    t1w_jsons = glob.glob(opj(session_dir, 'anat', '*_T1w.json'))
    for t1w_json in t1w_jsons:
        with open(t1w_json, 'r') as f:
            t1w_data = json.load(f)
        
        t1w_data['SkullStripped'] = False
        t1w_data['Manufacturer'] = 'Philips'
        t1w_data['MagneticFieldStrength'] = 7.0
        
        with open(t1w_json, 'w') as f:
            json.dump(t1w_data, f, indent=2)
    
    if t1w_jsons:
        print(f"  Updated {len(t1w_jsons)} T1w JSON files")
    print("✓ Anatomical metadata updated\n")

def main():
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <bids_dir> <subject> <session>")
        print("Example: python script.py /path/to/bids sub-01 ses-01")
        sys.exit(1)
    
    bids_dir = sys.argv[1]
    subject = sys.argv[2]
    session = sys.argv[3]
    
    # Construct session directory
    session_dir = opj(bids_dir, subject, session)
    
    # Check if session directory exists
    if not os.path.exists(session_dir):
        print(f"Error: Session directory does not exist: {session_dir}")
        sys.exit(1)
    
    print(f"Processing: {subject}/{session}")
    print(f"Session directory: {session_dir}\n")
    
    # Step 1: Update functional metadata
    print("Step 1: Updating functional JSON metadata...")
    update_functional_metadata(session_dir, subject)
    print("✓ Functional metadata updated\n")
    
    # Step 2: Fix TR in NIfTI headers
    print("Step 2: Fixing TR in NIfTI headers...")
    fix_nifti_tr(session_dir, subject)
    print("✓ TR fixed in NIfTI headers\n")
    
    # Step 3: Update EPI fieldmap metadata
    print("Step 3: Updating EPI fieldmap metadata...")
    update_epi_fieldmap_metadata(session_dir, subject)
    print("✓ EPI fieldmap metadata updated\n")
    
    # Step 4: Update legacy fieldmap metadata (if exists)
    print("Step 4: Updating legacy fieldmap metadata...")
    update_fieldmap_metadata(session_dir, subject)
    print("✓ Legacy fieldmap metadata updated\n")
    
    # Step 5: Update anatomical metadata
    print("Step 5: Updating anatomical metadata...")
    
    print("All steps completed successfully!")

if __name__ == "__main__":
    main()