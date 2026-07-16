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
import pandas as pd

def correct_task_columns_event(base_dir, subject, session):
    # Find events files
    event_files = glob.glob(
        os.path.join(base_dir, subject, session, "func", "*_events.tsv")
    )
    for event_file in event_files:
        # Load file
        event_table = pd.read_csv(event_file, sep="\t")
        event_table = event_table.loc[:, ~event_table.columns.str.contains("Unnamed")]

        # Rename column if it exists
        if 'task' in event_table.columns:
            event_table = event_table.rename(columns={'task': 'task_condition'})

        # Ensure integer-typed columns stay integer despite NaNs
        for col in ['stim_stair_val', 'response_val']:
            if col in event_table.columns:
                event_table[col] = event_table[col].astype('Int64')

        # Save back
        event_table.to_csv(event_file, sep="\t", index=False, na_rep='n/a')
        print(f"correction done for {event_file}")
        
def fix_skullstripped_key(json_dir):

    """

    Add SkullStripped only to relevant imaging JSON files

    and remove it from non-relevant ones inside a directory.

    """

    if not os.path.exists(json_dir):

        return

    keep_suffixes = (

        "_T1w.json",

        "_T2w.json",

        "_bold.json",

        "_sbref.json",
        
        "_magnitude.json",

    )

    for f in os.listdir(json_dir):

        if not f.endswith(".json"):

            continue

        json_path = os.path.join(json_dir, f)

        with open(json_path, "r") as jf:

            data = json.load(jf)

        is_valid = f.endswith(keep_suffixes)

        changed = False

        # add if needed

        if is_valid:

            if "SkullStripped" not in data:

                data["SkullStripped"] = False

                changed = True

                print(f"Added 'SkullStripped' key to {json_path}")

        # remove if not needed

        else:

            if "SkullStripped" in data:

                del data["SkullStripped"]

                changed = True

                print(f"Removed 'SkullStripped' key from {json_path}")

        if changed:

            with open(json_path, "w") as jf:

                json.dump(data, jf, indent=4)

def bidsify_func(sub, ses, base_dir):
    """
    Apply BIDS-related JSON fixes to functional data.
    """
    func_dir = os.path.join(base_dir, sub, ses, "func")
    fix_skullstripped_key(func_dir)

def bidsify_anat(sub, ses, base_dir):
    """
    Apply BIDS-related JSON fixes to anatomical data.
    """
    anat_dir = os.path.join(base_dir, sub, ses, "anat")
    fix_skullstripped_key(anat_dir)

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