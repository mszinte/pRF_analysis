#!/usr/bin/env python3
"""
Check integrity of NIfTI files (.nii / .nii.gz).

This script attempts to fully load each NIfTI file using nibabel.
If the file is corrupted or truncated, an exception will be raised.

Usage:
    python check_nifti_integrity.py /path/to/bids_or_sourcedata
"""

import sys
from pathlib import Path
import nibabel as nib


def check_nifti_file(nifti_path):
    """
    Try to fully load a NIfTI file.
    Raises an exception if the file is corrupted.
    """
    img = nib.load(str(nifti_path))
    _ = img.get_fdata()  # force full data read


def main(root_dir):
    root_dir = Path(root_dir)

    nifti_files = list(root_dir.rglob("*.nii")) + list(root_dir.rglob("*.nii.gz"))

    if not nifti_files:
        print("No NIfTI files found in {}".format(root_dir))
        return

    print("Found {} NIfTI files\n".format(len(nifti_files)))

    corrupted = []

    for nifti in sorted(nifti_files):
        try:
            check_nifti_file(nifti)
            print("[OK] {}".format(nifti))
        except Exception as exc:
            print("[CORRUPTED] {} -> {}".format(nifti, exc))
            corrupted.append(nifti)

    print("\n==============================")
    print("Summary")
    print("==============================")
    print("Total files checked : {}".format(len(nifti_files)))
    print("Corrupted files     : {}".format(len(corrupted)))

    if corrupted:
        print("\nList of corrupted files:")
        for f in corrupted:
            print(" - {}".format(f))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_nifti_integrity.py /path/to/data")
        sys.exit(1)

    main(sys.argv[1])