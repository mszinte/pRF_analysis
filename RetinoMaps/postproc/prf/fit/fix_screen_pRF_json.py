"""
-----------------------------------------------------------------------------------------
fix_screen_pRF_json.py
-----------------------------------------------------------------------------------------
Goal of the script:
[CORRECTION #29] Les events.json de la tache pRF contiennent la geometrie de l'ecran de la
salle d'eye-tracking (Display++ 0.696 x 0.391 m, avec 1.20 m) au lieu de celle de l'IRM
(ProPixx 0.77333 x 0.435 m a 1.20 m, pRFexp/config/scrConfig.m avec const.room = 1).
Le fit pRF (screen_utils.get_screen_settings) lit ces json : il a suppose 18.51 deg de
hauteur au lieu de 20.55 deg. Ce script remet la geometrie de l'IRM dans les json.
Pour chaque json, la geometrie est lue si possible dans le matFile Matlab du run
(sourcedata/<sub>/<ses>/add/<sub>_<ses>_task-<task>_run-0N_matFile.mat, config.scr :
disp_sizeX/Y en mm, dist en cm, scr_sizeX/Y en px) ; si le matFile indique la salle 2
(696 x 391 mm ou 108 cm) le json n'est PAS ecrit (ERREUR affichee). Sans matFile, les
constantes de l'IRM sont utilisees (0.77333 x 0.435 m, 1.20 m, 1920 x 1080 px).
Par defaut rien n'est ecrit (verification) ; --apply ecrit (copie <json>.bak_screen
avant ecriture, jamais ecrasee), refuse les liens symboliques (git-annex : faire
`datalad unlock` avant) puis verifie le resultat avec screen_utils.get_screen_settings.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
--apply: write the json files (default: dry-run, nothing written)
--all-tasks: also fix the events.json of every task run in the MRI
             (rest, pRF, PurLoc, SacLoc, PurVELoc, SacVELoc, pMF) for BIDS consistency
-----------------------------------------------------------------------------------------
Output(s):
corrected *_events.json (StimulusPresentation: ScreenSize, ScreenDistance, ScreenResolution)
backup copies *_events.json.bak_screen
summary table printed in the terminal
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/prf/fit
2. (datalad dataset) unlock the json files first
>> datalad unlock /scratch/mszinte/data/RetinoMaps/sub-*/ses-*/func/*_events.json
3. run python command
>> python fix_screen_pRF_json.py [main directory] [project name] [--apply] [--all-tasks]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/prf/fit
python fix_screen_pRF_json.py /scratch/mszinte/data RetinoMaps                        # verification (pRF)
python fix_screen_pRF_json.py /scratch/mszinte/data RetinoMaps --apply                # ecriture (pRF)
python fix_screen_pRF_json.py /scratch/mszinte/data RetinoMaps --all-tasks            # verification (toutes les taches IRM)
python fix_screen_pRF_json.py /scratch/mszinte/data RetinoMaps --all-tasks --apply    # ecriture (toutes les taches IRM)
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net) with Claude
-----------------------------------------------------------------------------------------
"""
# General imports
import os
import re
import sys
import glob
import json
import shutil
import argparse

# Inputs
parser = argparse.ArgumentParser(description="Fix screen geometry in BIDS events.json (MRI screen)")
parser.add_argument("main_dir", help="main project directory")
parser.add_argument("project_dir", help="project name (correspond to directory)")
parser.add_argument("--apply", action="store_true", help="write the json files (default: dry-run)")
parser.add_argument("--all-tasks", action="store_true", help="fix the events.json of all MRI tasks")
args = parser.parse_args()
main_dir = args.main_dir
project_dir = args.project_dir
root = f"{main_dir}/{project_dir}"

# MRI settings (pRFexp/config/scrConfig.m, const.room = 1 ; expLauncher.m const.desiredRes)
screen_size_m = [0.77333, 0.435]
screen_distance_m = 1.20
screen_res = [1920, 1080]

# Eye-tracking room settings (scrConfig.m, const.room = 2) : must never end up in the json
room2_size_mm = [696, 391]
room2_dist_cm = 108

# Tasks
mri_tasks = ["rest", "pRF", "PurLoc", "SacLoc", "PurVELoc", "SacVELoc", "pMF"]
tasks = mri_tasks if args.all_tasks else ["pRF"]

# Find json files
json_fns = []
for task in tasks:
    # only BIDS *_events.json (the eye-tracking *_physioevents.json are left untouched)
    json_fns.extend(sorted(fn for fn in glob.glob(f"{root}/sub-*/**/*task-{task}_*events.json", recursive=True)
                           if os.path.basename(fn).endswith('_events.json')))
    root_json = f"{root}/task-{task}_events.json"
    if os.path.lexists(root_json):
        json_fns.append(root_json)
json_fns = list(dict.fromkeys(json_fns))  # remove duplicates, keep order
if not json_fns:
    sys.exit(f"ERROR: no events.json found for tasks {tasks} in {root}")


def find_matfile(json_fn):
    """
    Find the Matlab matFile matching a BIDS events.json (same subject, session, task, run).
    Returns (matfile path or None, entities dict).
    """
    base = os.path.basename(json_fn)
    ents = {k: (m.group(1) if m else None) for k, m in (
        ('sub', re.search(r'sub-([a-zA-Z0-9]+)', base)),
        ('ses', re.search(r'ses-([a-zA-Z0-9]+)', base)),
        ('task', re.search(r'task-([a-zA-Z0-9]+)', base)),
        ('run', re.search(r'run-0*(\d+)', base)))}
    if ents['sub'] is None:
        return None, ents
    sub, ses, task = f"sub-{ents['sub']}", f"ses-{ents['ses']}", ents['task']
    run_txt = f"run-{int(ents['run']):02d}" if ents['run'] else "run-*"
    # naming of pRFexp/config/dirSaveFile.m : data/<sub>/<ses>/add/<sub>_<ses>_task-<task>_run-0N_matFile.mat
    candidates = []
    if ents['ses']:
        candidates = sorted(glob.glob(f"{root}/sourcedata/{sub}/{ses}/add/{sub}_{ses}_task-{task}_{run_txt}_matFile.mat"))
    if not candidates:
        candidates = sorted(glob.glob(f"{root}/sourcedata/{sub}/**/*task-{task}*{run_txt}*matFile*.mat", recursive=True))
        if ents['ses']:
            same_ses = [c for c in candidates if f"_{ses}_" in os.path.basename(c)]
            candidates = same_ses if same_ses else candidates
    return (candidates[0] if candidates else None), ents


def read_matfile_screen(matfile):
    """
    Read config.scr of a Matlab matFile (as in postproc/pmf/fit/get_screen_info_mat.py).
    Returns size (mm), distance (cm) and resolution (px).
    """
    import scipy.io
    mat = scipy.io.loadmat(matfile)
    scr = mat['config'][0][0]['scr'][0][0]
    size_mm = [float(scr['disp_sizeX'][0][0]), float(scr['disp_sizeY'][0][0])]
    dist_cm = float(scr['dist'][0][0])
    res_px = [int(scr['scr_sizeX'][0][0]), int(scr['scr_sizeY'][0][0])]
    return size_mm, dist_cm, res_px


# Process json files
summary = []
n_error = 0
for json_fn in json_fns:
    rel_fn = os.path.relpath(json_fn, root)
    link_txt = ' [symlink: datalad unlock]' if os.path.islink(json_fn) else ''

    # New values : matFile if found, MRI constants otherwise
    matfile, ents = find_matfile(json_fn)
    source = 'constants'
    new_size, new_dist, new_res = screen_size_m, screen_distance_m, screen_res
    if matfile is not None:
        try:
            size_mm, dist_cm, res_px = read_matfile_screen(matfile)
        except Exception as e:
            print(f"WARNING: cannot read {matfile} ({e}) : MRI constants used")
        else:
            is_room2 = (abs(size_mm[0] - room2_size_mm[0]) < 1 and abs(size_mm[1] - room2_size_mm[1]) < 1) \
                       or abs(dist_cm - room2_dist_cm) < 0.5
            if is_room2:
                print(f"ERROR: {os.path.relpath(matfile, root)} indicates the eye-tracking room "
                      f"({size_mm[0]:g} x {size_mm[1]:g} mm, {dist_cm:g} cm) : {rel_fn} NOT written, check this run")
                summary.append((rel_fn, '-', '-', 'ERROR room 2 (matFile)'))
                n_error += 1
                continue
            source = 'matFile'
            new_size = [round(size_mm[0] / 1000, 5), round(size_mm[1] / 1000, 5)]
            new_dist = round(dist_cm / 100, 4)
            new_res = res_px
            if abs(new_size[0] - screen_size_m[0]) > 1e-3 or abs(new_size[1] - screen_size_m[1]) > 1e-3 \
               or abs(new_dist - screen_distance_m) > 1e-3 or list(new_res) != screen_res:
                print(f"WARNING: {os.path.relpath(matfile, root)} differs from MRI constants "
                      f"({new_size}, {new_dist}, {new_res}) : matFile values used")

    # Read json
    with open(json_fn, 'r', encoding='utf-8') as f:
        data = json.load(f)
    stim = data.setdefault('StimulusPresentation', {})
    old = (stim.get('ScreenSize'), stim.get('ScreenDistance'), stim.get('ScreenResolution'))
    new = (list(new_size), new_dist, list(new_res))
    old_txt = f"{old[0]}, {old[1]}, {old[2]}"
    new_txt = f"{new[0]}, {new[1]}, {new[2]}"
    if old == new:
        summary.append((rel_fn, old_txt, 'unchanged', source))
        continue
    stim['ScreenSize'] = new[0]
    stim['ScreenDistance'] = new[1]
    stim['ScreenResolution'] = new[2]

    # Write json
    if args.apply:
        if os.path.islink(json_fn):
            print(f"ERROR: {rel_fn} is a symbolic link (git-annex) : not written. "
                  f"Run `datalad unlock {json_fn}` then rerun with --apply")
            summary.append((rel_fn, old_txt, '-', 'ERROR symlink'))
            n_error += 1
            continue
        backup_fn = f"{json_fn}.bak_screen"
        if not os.path.exists(backup_fn):
            shutil.copy2(json_fn, backup_fn)
        with open(json_fn, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4, ensure_ascii=False)
    summary.append((rel_fn, old_txt, new_txt, source + link_txt))

# Summary table
w_file = max(len('file'), max(len(s[0]) for s in summary))
w_old = max(len('old values'), max(len(s[1]) for s in summary))
w_new = max(len('new values'), max(len(s[2]) for s in summary))
print(f"\n{'file':<{w_file}} | {'old values (size m, dist m, res px)':<{w_old}} | {'new values':<{w_new}} | source")
print(f"{'-' * w_file}-+-{'-' * w_old}-+-{'-' * w_new}-+-{'-' * 22}")
for rel_fn, old_txt, new_txt, source in summary:
    print(f"{rel_fn:<{w_file}} | {old_txt:<{w_old}} | {new_txt:<{w_new}} | {source}")
print(f"\n{len(summary)} json file(s), {n_error} error(s)")
print('WRITTEN (--apply)' if args.apply else 'Dry-run only : rerun with --apply to write')

# Verification with the function used by the pRF fits
if args.apply:
    # direct check of every json (get_screen_settings only reads the first json of each subject)
    print("\nDirect check of every json:")
    n_wrong = 0
    for json_fn in json_fns:
        with open(json_fn, 'r', encoding='utf-8') as f:
            stim = json.load(f).get('StimulusPresentation', {})
        size, dist = stim.get('ScreenSize', [0, 0]), stim.get('ScreenDistance', 0)
        if abs(size[0] - screen_size_m[0]) > 1e-3 or abs(size[1] - screen_size_m[1]) > 1e-3 \
           or abs(dist - screen_distance_m) > 1e-3:
            print(f"WRONG: {os.path.relpath(json_fn, root)}: {size}, {dist}")
            n_wrong += 1
    print(f"{len(json_fns) - n_wrong}/{len(json_fns)} json file(s) with the MRI screen geometry")

    sys.path.append(os.path.abspath("{}/../../../../analysis_code/utils".format(os.path.dirname(os.path.abspath(__file__)))))
    try:
        from screen_utils import get_screen_settings
        import bids  # noqa: F401 (pybids, used by get_screen_settings)
    except ImportError as e:
        print(f"WARNING: verification skipped, screen_utils/pybids not available ({e})")
    else:
        print("\nVerification with screen_utils.get_screen_settings (task pRF):")
        subjects = sorted({f"sub-{m.group(1)}" for fn in json_fns if '_task-pRF_' in os.path.basename(fn)
                           for m in [re.search(r'sub-([a-zA-Z0-9]+)', os.path.basename(fn))] if m})
        n_bad, n_unverified = 0, 0
        for subject in subjects:
            try:
                screen_size_cm, screen_distance_cm = get_screen_settings(main_dir, project_dir, subject[4:], 'pRF')
            except Exception as e:
                print(f"WARNING: {subject}: NOT VERIFIED, get_screen_settings failed ({e})")
                n_unverified += 1
                continue
            ok = abs(screen_size_cm[0] - 77.33) < 0.01 and abs(screen_size_cm[1] - 43.5) < 0.01 \
                 and abs(screen_distance_cm - 120) < 0.01
            print(f"{subject}: screen size {screen_size_cm} cm, distance {screen_distance_cm} cm -> {'OK' if ok else 'WRONG'}")
            n_bad += not ok
        assert n_bad == 0, f"{n_bad} subject(s) without the MRI screen geometry (77.33 x 43.5 cm at 120 cm)"
        if n_unverified:
            print(f"WARNING: {n_unverified} subject(s) not verified (see above)")
        else:
            print("All subjects: 77.33 x 43.5 cm at 120 cm")
