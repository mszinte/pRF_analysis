"""
-----------------------------------------------------------------------------------------
bids_utils.py
-----------------------------------------------------------------------------------------
Goal of the script:
Utility functions for BIDS conversion of amblyopia 7T pRF data
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Claude (Anthropic)
-----------------------------------------------------------------------------------------
"""

import os
import re
import json
import glob
import subprocess
import nibabel as nib
import numpy as np
opj = os.path.join


def maxk(A, k):
    """Get k largest values from array"""
    B = np.sort(np.ndarray.flatten(A))
    t = B.size
    y = B[t-k:]
    return y


def generate_psir_images(dcm2niix_dir, prefix="PARREC_WIPPSIR1mm3SENSE"):
    """
    Generate PSIR reconstructed images from PARREC converted files.
    Creates 4 outputs: MP2RAGE_Mougin, PSIR_Mougin, PSIRbeta_Mougin, TI2_PSIR_Mougin
    
    Based on Mougin et al. (2016) MRM 76:1512-1516
    
    Parameters
    ----------
    dcm2niix_dir : str
        Directory with dcm2niix output containing PSIR files
    prefix : str
        Prefix pattern for PSIR files (default: "PARREC_WIPPSIR1mm3SENSE")
        
    Returns
    -------
    output_files : dict
        Dictionary with paths to the 4 generated files
    """
    # Find PSIR files with the prefix pattern
    psir_files = glob.glob(opj(dcm2niix_dir, f"{prefix}*_t787.nii.gz"))
    if not psir_files:
        raise FileNotFoundError(f"No PSIR files found with prefix {prefix} in {dcm2niix_dir}")
    
    # Extract the actual prefix from the first file
    actual_prefix = os.path.basename(psir_files[0]).replace('_t787.nii.gz', '')
    
    print(f"  Loading PSIR data with prefix: {actual_prefix}")
    
    # Load magnitude images
    imgm1 = nib.load(opj(dcm2niix_dir, f"{actual_prefix}_t787.nii.gz"))
    imgm2 = nib.load(opj(dcm2niix_dir, f"{actual_prefix}_t2387.nii.gz"))
    
    # Load phase images
    imgp1 = nib.load(opj(dcm2niix_dir, f"{actual_prefix}_ph_t787.nii.gz"))
    imgp2 = nib.load(opj(dcm2niix_dir, f"{actual_prefix}_ph_t2387.nii.gz"))
    
    # Get data
    imgm1_data = np.squeeze(imgm1.get_fdata())
    imgm2_data = np.squeeze(imgm2.get_fdata())
    imgp1_data = np.squeeze(imgp1.get_fdata())
    imgp2_data = np.squeeze(imgp2.get_fdata())
    
    print(f"  Data shape: {imgm1_data.shape}")
    
    # Step 1: Create complex images: S = magnitude × exp(i×phase)
    S1 = np.multiply(imgm1_data, np.exp(1j * imgp1_data))
    S2 = np.multiply(imgm2_data, np.exp(1j * imgp2_data))
    
    # Step 2: Compute polarity from phase
    T1w = np.real(np.multiply(np.conj(S1), S2))
    f = np.sign(T1w)  # Polarity map: +1 or -1
    
    # Apply polarity to magnitude of S1
    S1m = np.multiply(np.abs(S1), f)
    
    # Normalization image (sum of magnitudes)
    S0 = np.abs(S1m) + np.abs(S2)
    
    # Step 3: Compute MP2RAGE
    MP2RAGEn = np.real(np.multiply(np.conj(S1), S2))
    MP2RAGEd = np.square(np.abs(S1)) + np.square(np.abs(S2))
    MP2RAGE = np.divide(MP2RAGEn, MP2RAGEd, out=np.zeros_like(MP2RAGEn), where=(MP2RAGEd != 0))
    
    # Step 4: Compute beta parameters (regularization)
    aa = 1
    beta = aa * np.min(maxk(np.abs(S2), round(S2.size / 100)))
    
    # Step 5: Compute PSIRbeta (with beta regularization, scaled to 0-400 range)
    PSIRbeta = (np.divide(S1m - beta, S0 + beta, out=np.zeros_like(S1m), where=(S0 + beta != 0)) + 1) * 200
    
    # Step 6: Compute simple PSIR (without beta)
    PSIR = np.divide(S1m, S0, out=np.zeros_like(S1m), where=(S0 != 0))
    
    # Step 7: Compute TI2_PSIR (PSIR scaled by magnitude of S2)
    TI2_PSIR = np.multiply(PSIR + 1, np.abs(S2))
    
    # Save all outputs
    output_files = {}
    
    # Save MP2RAGE
    MP2RAGE_img = nib.Nifti1Image(MP2RAGE, imgm1.affine, imgm1.header)
    MP2RAGE_img.set_data_dtype(np.float32)
    MP2RAGE_output = opj(dcm2niix_dir, f"{actual_prefix}_MP2RAGE_Mougin.nii.gz")
    nib.save(MP2RAGE_img, MP2RAGE_output)
    output_files['MP2RAGE'] = MP2RAGE_output
    print(f"  Saved: {os.path.basename(MP2RAGE_output)}")
    
    # Save PSIR
    PSIR_img = nib.Nifti1Image(PSIR, imgm1.affine, imgm1.header)
    PSIR_img.set_data_dtype(np.float32)
    PSIR_output = opj(dcm2niix_dir, f"{actual_prefix}_PSIR_Mougin.nii.gz")
    nib.save(PSIR_img, PSIR_output)
    output_files['PSIR'] = PSIR_output
    print(f"  Saved: {os.path.basename(PSIR_output)}")
    
    # Save PSIRbeta (main output for T1w)
    PSIRbeta_img = nib.Nifti1Image(PSIRbeta, imgm1.affine, imgm1.header)
    PSIRbeta_img.set_data_dtype(np.float32)
    PSIRbeta_output = opj(dcm2niix_dir, f"{actual_prefix}_PSIRbeta_Mougin.nii.gz")
    nib.save(PSIRbeta_img, PSIRbeta_output)
    output_files['PSIRbeta'] = PSIRbeta_output
    print(f"  Saved: {os.path.basename(PSIRbeta_output)}")
    
    # Save TI2_PSIR
    TI2_PSIR_img = nib.Nifti1Image(TI2_PSIR, imgm1.affine, imgm1.header)
    TI2_PSIR_img.set_data_dtype(np.float32)
    TI2_PSIR_output = opj(dcm2niix_dir, f"{actual_prefix}_TI2_PSIR_Mougin.nii.gz")
    nib.save(TI2_PSIR_img, TI2_PSIR_output)
    output_files['TI2_PSIR'] = TI2_PSIR_output
    print(f"  Saved: {os.path.basename(TI2_PSIR_output)}")
    
    return output_files


def parse_readme(readme_path):
    """
    Parse Readme.txt to extract scan information and create task mappings.
    
    Parameters
    ----------
    readme_path : str
        Path to Readme.txt file
        
    Returns
    -------
    scans : list of dict
        List of scan information with task mappings
    """
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
    
    return scans


def create_verification_file(scans, subject, readme_path, output_path):
    """
    Create human-readable verification file for BIDS mapping.
    
    Parameters
    ----------
    scans : list of dict
        List of scan information
    subject : str
        Subject BIDS code
    readme_path : str
        Path to original Readme.txt
    output_path : str
        Path to save verification file
    """
    with open(output_path, 'w') as f:
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


def run_dcm2niix(parrec_dir, output_dir):
    """
    Run dcm2niix to convert PARREC files to NIfTI.
    Skips conversion if output directory exists and is not empty.
    
    Parameters
    ----------
    parrec_dir : str
        Directory containing PARREC files
    output_dir : str
        Output directory for NIfTI files
        
    Returns
    -------
    bool
        True if conversion was run, False if skipped
    """
    # Check if output directory exists and is not empty
    if os.path.exists(output_dir) and os.listdir(output_dir):
        print(f"  dcm2niix output directory already exists and is not empty: {output_dir}")
        print(f"  Skipping dcm2niix conversion")
        return False
    
    # Create output directory and run conversion
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"dcm2niix -o {output_dir} -z y {parrec_dir}"
    subprocess.run(cmd, shell=True, check=True)
    return True


def copy_anatomical(dcm2niix_dir, session_dir, subject):
    """
    Copy anatomical T1w files to BIDS structure.
    
    Parameters
    ----------
    dcm2niix_dir : str
        Directory with dcm2niix output
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    """
    # Find T1w files (PSIRbeta_Mougin output)
    t1w_nii = glob.glob(opj(dcm2niix_dir, '*PSIRbeta_Mougin.nii.gz'))
    t1w_json = glob.glob(opj(dcm2niix_dir, '*PSIR*.json'))
    
    if t1w_nii:
        # Copy NIfTI
        dst_nii = opj(session_dir, 'anat', f'{subject}_ses-01_T1w.nii.gz')
        subprocess.run(f"cp {t1w_nii[0]} {dst_nii}", shell=True, check=True)
        
        # Copy JSON
        if t1w_json:
            dst_json = opj(session_dir, 'anat', f'{subject}_ses-01_T1w.json')
            subprocess.run(f"cp {t1w_json[0]} {dst_json}", shell=True, check=True)
        
        print(f"  T1w: {os.path.basename(t1w_nii[0])} -> {subject}_ses-01_T1w.nii.gz")
    else:
        print("  WARNING: T1w file not found")


def copy_functional_runs(scans, dcm2niix_dir, session_dir, subject, copy_commands_file):
    """
    Copy functional runs to BIDS structure based on scan mapping.
    
    Parameters
    ----------
    scans : list of dict
        List of scan information
    dcm2niix_dir : str
        Directory with dcm2niix output
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    copy_commands_file : str
        Path to save copy commands script
    """
    # Find all dcm2niix functional files
    dcm2niix_files = {}
    for filename in os.listdir(dcm2niix_dir):
        if filename.endswith('.nii.gz') and 'fMRI' in filename:
            # Extract the file number from the end of filename
            parts = filename.replace('.nii.gz', '').split('_')
            if parts[-1].isdigit():
                file_num = parts[-1]
                dcm2niix_files[file_num] = filename.replace('.nii.gz', '')
    
    print("  Found dcm2niix files:")
    for num, base in sorted(dcm2niix_files.items()):
        print(f"    File #{num}: {base}")
    
    # Create copy commands
    copy_commands = []
    copy_commands.append("#!/bin/bash\n")
    
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
            src_nii = opj(dcm2niix_dir, f"{source_base}.nii.gz")
            dst_nii = opj(session_dir, 'func', f"{target_base}.nii.gz")
            copy_commands.append(f"cp {src_nii} {dst_nii}\n")
            
            # JSON file
            src_json = opj(dcm2niix_dir, f"{source_base}.json")
            dst_json = opj(session_dir, 'func', f"{target_base}.json")
            copy_commands.append(f"cp {src_json} {dst_json}\n")
            
            print(f"  Scan {scan['scan_num']}: {source_base} -> {target_base}")
        else:
            print(f"  WARNING: File #{file_num} not found for scan {scan['scan_num']}")
    
    # Write and execute copy commands
    with open(copy_commands_file, 'w') as f:
        f.writelines(copy_commands)
    
    os.chmod(copy_commands_file, 0o755)
    subprocess.run(copy_commands_file, shell=True, check=True)


def copy_fieldmaps(dcm2niix_dir, session_dir, subject):
    """
    Copy fieldmap files to BIDS structure.
    
    Parameters
    ----------
    dcm2niix_dir : str
        Directory with dcm2niix output
    session_dir : str
        Session directory in BIDS structure
    subject : str
        Subject BIDS code
    """
    # Find fieldmap files
    fmap_phase_nii = glob.glob(opj(dcm2niix_dir, '*t2s_phase*_ph.nii.gz'))
    fmap_phase_json = glob.glob(opj(dcm2niix_dir, '*t2s_phase*_ph.json'))
    fmap_mag_nii = glob.glob(opj(dcm2niix_dir, '*t2s_phase*.nii.gz'))
    fmap_mag_nii = [f for f in fmap_mag_nii if not f.endswith('_ph.nii.gz')]
    
    if fmap_phase_nii:
        dst = opj(session_dir, 'fmap', f'{subject}_ses-01_fieldmap.nii.gz')
        subprocess.run(f"cp {fmap_phase_nii[0]} {dst}", shell=True)
        print(f"  Fieldmap phase: {os.path.basename(fmap_phase_nii[0])}")
    
    if fmap_phase_json:
        dst = opj(session_dir, 'fmap', f'{subject}_ses-01_fieldmap.json')
        subprocess.run(f"cp {fmap_phase_json[0]} {dst}", shell=True)
    
    if fmap_mag_nii:
        dst = opj(session_dir, 'fmap', f'{subject}_ses-01_magnitude.nii.gz')
        subprocess.run(f"cp {fmap_mag_nii[0]} {dst}", shell=True)
        print(f"  Fieldmap magnitude: {os.path.basename(fmap_mag_nii[0])}")


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
        data['RepetitionTime'] = 2.0
        data['SkullStripped'] = True
        data['PhaseEncodingDirection'] = 'i-'
        data['Manufacturer'] = 'Philips'
        data['ManufacturersModelName'] = 'Intera Achieva'
        data['MagneticFieldStrength'] = 7.0
        
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
        zooms[3] = 2.0
        hdr.set_zooms(zooms)
        new_img = nib.Nifti1Image(img.get_fdata(), img.affine, hdr)
        nib.save(new_img, func_file)
    
    print(f"  Fixed TR in {len(func_files)} NIfTI files")


def update_fieldmap_metadata(session_dir, subject):
    """
    Update fieldmap JSON with required BIDS metadata and IntendedFor field.
    
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
    
    # Update fieldmap JSON
    fieldmap_json = opj(session_dir, 'fmap', f'{subject}_ses-01_fieldmap.json')
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
        
        print("  Updated fieldmap JSON")
    
    # Update T1w JSON
    t1w_json = opj(session_dir, 'anat', f'{subject}_ses-01_T1w.json')
    if os.path.exists(t1w_json):
        with open(t1w_json, 'r') as f:
            t1w_data = json.load(f)
        
        t1w_data['SkullStripped'] = False
        t1w_data['Manufacturer'] = 'Philips'
        t1w_data['ManufacturersModelName'] = 'Intera Achieva'
        t1w_data['MagneticFieldStrength'] = 7.0
        
        with open(t1w_json, 'w') as f:
            json.dump(t1w_data, f, indent=2)
        
        print("  Updated T1w JSON")


def update_dataset_description(base_path):
    """
    Update dataset_description.json with required BIDS metadata.
    
    Parameters
    ----------
    base_path : str
        Base path of the dataset
    """
    dataset_json = opj(base_path, 'dataset_description.json')
    
    try:
        with open(dataset_json, 'r') as f:
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
    
    with open(dataset_json, 'w') as f:
        json.dump(data, f, indent=2)
