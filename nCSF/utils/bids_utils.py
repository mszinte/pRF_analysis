"""
-----------------------------------------------------------------------------------------
bids_utils.py
-----------------------------------------------------------------------------------------
Goal of the script:
Utility functions for BIDS conversion of nCSF data
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

import os
import json
import subprocess


def fast_copy(src, dst):
    """
    Copy a file or directory quickly using rsync with archive, compression, update, verbose.
    Handles both files and directories.
    """
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source not found: {src}")

    dst_dir = os.path.dirname(dst) if os.path.isfile(src) else dst
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    # Remove trailing slash for files, keep for directories
    if os.path.isfile(src):
        cmd = ["rsync", "-azuv", src, dst]
    else:
        cmd = ["rsync", "-azuv", src + "/", dst + "/"]

    subprocess.run(cmd, check=True)

def bidsify_fmap(sub, ses, base_dir):
    fmap_dir = os.path.join(base_dir, sub, ses, "fmap")
    func_dir = os.path.join(base_dir, sub, ses, "func")

    # 1. Rename TB1TFL → TB1map
    for f in os.listdir(fmap_dir):
        if "TB1TFL" in f:
            old_path = os.path.join(fmap_dir, f)
            new_name = f.replace("TB1TFL", "TB1map")
            new_path = os.path.join(fmap_dir, new_name)
            os.rename(old_path, new_path)

    # 2. Rename AP/PA fieldmaps robustly
    for f in os.listdir(fmap_dir):
        if "epi" in f and ("AP" in f or "PA" in f):
            old_path = os.path.join(fmap_dir, f)

            if "AP" in f:
                dirtag = "AP"
            else:
                dirtag = "PA"

            if f.endswith(".nii.gz"):
                ext = ".nii.gz"
            else:
                ext = ".json"

            new_name = "{}_{}_dir-{}_epi{}".format(sub, ses, dirtag, ext)
            new_path = os.path.join(fmap_dir, new_name)
            os.rename(old_path, new_path)

    # 3. Build IntendedFor list automatically from func/
    intended_for = []
    for f in sorted(os.listdir(func_dir)):
        if f.endswith("_bold.nii.gz"):
            intended_for.append("{}/func/{}".format(ses, f))

    # 4. Update JSON AP and PA
    for dir_tag in ["AP", "PA"]:
        json_path = os.path.join(
            fmap_dir,
            "{}_{}_dir-{}_epi.json".format(sub, ses, dir_tag)
        )

        if os.path.exists(json_path):
            with open(json_path, "r") as f:
                data = json.load(f)

            data["IntendedFor"] = intended_for

            with open(json_path, "w") as f:
                json.dump(data, f, indent=2)

    print("{} {} fmap BIDS-ified successfully.".format(sub, ses))

def bidsify_anat(sub, ses, base_dir):
    anat_dir = os.path.join(base_dir, sub, ses, "anat")

    # Identify files
    files = os.listdir(anat_dir)

    # Remove forbidden diffusion files
    for f in files:
        if f.endswith(".bval") or f.endswith(".bvec"):
            os.remove(os.path.join(anat_dir, f))

    # Process MP2RAGE components
    for f in files:
        old = os.path.join(anat_dir, f)

        # Skip removed files
        if not os.path.exists(old):
            continue

        ext = ".nii.gz" if f.endswith(".nii.gz") else ".json"

        # UNIDEN -> T1w
        if "UNIDEN" in f:
            new = "{}_{}_T1w{}".format(sub, ses, ext)
            os.rename(old, os.path.join(anat_dir, new))
            continue

        # UNI (non denoised)
        if "UNI" in f and "UNIDEN" not in f:
            new = "{}_{}_acq-UNI_MP2RAGE{}".format(sub, ses, ext)
            os.rename(old, os.path.join(anat_dir, new))
            continue

        # Inversion 1
        if "inv-1" in f and "part-mag" in f:
            new = "{}_{}_inv-1_part-mag_MP2RAGE{}".format(sub, ses, ext)
            os.rename(old, os.path.join(anat_dir, new))
            continue

        # Inversion 2
        if "inv-2" in f and "part-mag" in f:
            new = "{}_{}_inv-2_part-mag_MP2RAGE{}".format(sub, ses, ext)
            os.rename(old, os.path.join(anat_dir, new))
            continue

    print("Anat for {} {} BIDS-ified successfully.".format(sub, ses))