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
# Debug
import ipdb
deb = ipdb.set_trace

import os
import glob
import json
import scipy.io
import pandas as pd


def clean_anat_mp2rage(anat_dir):
    """
    Keep only MP2RAGE UNIDEN T1w image and remove all other MP2RAGE-related files.

    This function assumes that original files are safely stored in sourcedata.
    """
    if not os.path.isdir(anat_dir):
        return

    for fname in os.listdir(anat_dir):

        # Work only on MP2RAGE-related files
        if "MP2RAGE" not in fname and "UNIDEN" not in fname:
            continue

        fpath = os.path.join(anat_dir, fname)

        # Keep UNIDEN T1w (nii.gz and json)
        if (
            "UNIDEN" in fname
            and "_T1w" in fname
            and (fname.endswith(".nii.gz") or fname.endswith(".json"))
        ):
            continue

        # Remove everything else related to MP2RAGE
        if os.path.isfile(fpath):
            os.remove(fpath)

    print("Kept only UNIDEN T1w and removed other MP2RAGE files in {}".format(anat_dir))


def update_mp2rage_json(anat_dir, tr_exc=0, tr_prep=5, nshots=1):
    """
    Update MP2RAGE JSON files with required BIDS fields.
    tr_exc is intentionally left to 0 (temporary placeholder).
    """
    if not os.path.exists(anat_dir):
        return

    for f in os.listdir(anat_dir):
        if f.endswith("MP2RAGE.nii.gz"):
            json_path = os.path.join(anat_dir, f.replace(".nii.gz", ".json"))

            if os.path.exists(json_path):
                with open(json_path, "r") as jf:
                    data = json.load(jf)
            else:
                data = {}

            data["RepetitionTimeExcitation"] = tr_exc
            data["RepetitionTimePreparation"] = tr_prep
            data["NumberShots"] = nshots
            data["SkullStripped"] = False

            with open(json_path, "w") as jf:
                json.dump(data, jf, indent=4)

            print("Updated MP2RAGE JSON for {}".format(f))

def correct_task_columns_event(base_dir, subject, session):
    # Find DeepMReye event files
    event_files = glob.glob(os.path.join(base_dir, subject, session, "func", '*task-DeepMReye*_events.tsv'))

    for event_file in event_files:
        # Load the file
        event_table = pd.read_csv(event_file, sep="\t", keep_default_na=False)
        
        # Rename the 'task' column to 'trial_type' if it exists
        if 'task' in event_table.columns:
            event_table.rename(columns={'task': 'trial_type'}, inplace=True)

        # Convert columns that should be integers to nullable Int64
        for col in event_table.columns:
            if event_table[col].dtype == float and event_table[col].dropna().apply(float.is_integer).all():
                event_table[col] = event_table[col].astype('Int64')

        # Save the file
        event_table.to_csv(event_file, sep="\t", index=False, na_rep='n/a')
        print("correction done for task column in {}".format(event_file))
        
def create_subject_task_events_json(base_dir, subject, session, task, run="01"):
    """
    Create the *_events.json sidecar associated with an existing *_events.tsv file.
    """
    source_func_dir = os.path.join(base_dir, "sourcedata", subject, session, "func")
    task_json_dir = os.path.join(base_dir, subject, session)


    mat_file_fn = os.path.join(
        source_func_dir,
        "{}_{}_task-{}_dir-PA_run-{}_matFile.mat".format(subject, session, task, run)
    )

    if not os.path.exists(mat_file_fn):
        print("WARNING: missing matlab file {}".format(mat_file_fn))
        return

    mat = scipy.io.loadmat(mat_file_fn)
    config = mat["config"][0, 0]
    scr = config["scr"][0, 0]

    task_metadata = {
        "TaskName": task,
        "StimulusPresentation": {
            "ScreenDistance": float(scr["dist"][0][0] / 100),
            "ScreenOrigin": ["top", "left"],
            "ScreenRefreshRate": int(scr["hz"][0][0]),
            "ScreenSize": [
                float(scr["disp_sizeX"][0][0] / 1000),
                float(scr["disp_sizeY"][0][0] / 1000)
            ],
            "ScreenResolution": [
                int(scr["scr_sizeX"][0][0]),
                int(scr["scr_sizeY"][0][0])
            ]
        },
        "SoftwareName": "Psychtoolbox",
        "SoftwareVersion": ""
    }

    task_json_fn = os.path.join(
        task_json_dir,
        "{}_{}_task-{}_events.json".format(task, subject, session)
    )

    with open(task_json_fn, "w") as f:
        json.dump(task_metadata, f, indent=4)

    print("Created events JSON for {} {} task {} run {}".format(
        subject, session, task, run)
    )


def add_keys_json_func(json_dir):
    """
    Add missing Instructions and TaskDescription fields to BOLD JSON files.
    """
    if not os.path.exists(json_dir):
        return

    for f in os.listdir(json_dir):
        if f.endswith(".json") and "bold" in f:
            json_path = os.path.join(json_dir, f)

            with open(json_path, "r") as jf:
                data = json.load(jf)

            if "Instructions" not in data:
                data["Instructions"] = (
                    "keep fixation on the center and report noise orientation "
                    "by clicking with left or right thumb"
                )

            if "TaskDescription" not in data:
                data["TaskDescription"] = (
                    "keep fixation on the center and report noise orientation "
                    "by clicking with left or right thumb"
                )

            with open(json_path, "w") as jf:
                json.dump(data, jf, indent=4)


def add_skullstripped_key(json_dir):
    """
    Add SkullStripped field to JSON files if missing.
    """
    if not os.path.exists(json_dir):
        return

    for f in os.listdir(json_dir):
        if f.endswith(".json"):
            json_path = os.path.join(json_dir, f)

            with open(json_path, "r") as jf:
                data = json.load(jf)

            if "SkullStripped" not in data:
                data["SkullStripped"] = False
                with open(json_path, "w") as jf:
                    json.dump(data, jf, indent=4)

                print("Added 'SkullStripped' key to {}".format(json_path))


def bidsify_func(sub, ses, base_dir):
    """
    Apply BIDS-related JSON fixes to functional data.
    """
    func_dir = os.path.join(base_dir, sub, ses, "func")
    add_skullstripped_key(func_dir)
    add_keys_json_func(func_dir)


def bidsify_anat(sub, ses, base_dir):
    """
    Apply BIDS-related JSON fixes to anatomical data.
    """
    anat_dir = os.path.join(base_dir, sub, ses, "anat")
    clean_anat_mp2rage(anat_dir)
    update_mp2rage_json(anat_dir)
    add_skullstripped_key(anat_dir)

    print("Anat BIDS-ified for {} {}".format(sub, ses))


def bidsify_fmap(sub, ses, base_dir):
    """
    Add IntendedFor and SkullStripped fields to fmap JSON files.
    """
    fmap_dir = os.path.join(base_dir, sub, ses, "fmap")
    func_dir = os.path.join(base_dir, sub, ses, "func")

    if not os.path.exists(fmap_dir) or not os.path.exists(func_dir):
        return

    intended_for = [
        "{}/func/{}".format(ses, f)
        for f in sorted(os.listdir(func_dir))
        if f.endswith("_bold.nii.gz")
    ]

    for f in os.listdir(fmap_dir):
        if f.endswith(".json"):
            json_path = os.path.join(fmap_dir, f)

            with open(json_path, "r") as jf:
                data = json.load(jf)

            data["IntendedFor"] = intended_for

            if "SkullStripped" not in data:
                data["SkullStripped"] = False

            with open(json_path, "w") as jf:
                json.dump(data, jf, indent=4)

    print("Fmap BIDS-ified for {} {}".format(sub, ses))