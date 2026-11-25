#!/usr/bin/env python3
"""Rename every <pp_data>/<subject or group>/<format>/prfmodel folder to prfmodel_old.

Usage:
    python rename_prfmodel.py <pp_data_dir>          # dry run: only prints what would be done
    python rename_prfmodel.py <pp_data_dir> --run    # actually renames
Example (cluster):
    python rename_prfmodel.py /scratch/mszinte/data/amblyo7T_prf/derivatives/pp_data
    python rename_prfmodel.py /scratch/mszinte/data/amblyo7T_prf/derivatives/pp_data --run
"""
import os
import sys
import glob

if len(sys.argv) < 2:
    sys.exit(__doc__)

pp_dir = sys.argv[1]
do_run = '--run' in sys.argv
old_name, new_name = 'prfmodel', 'prfmodel_old'

targets = sorted(glob.glob(os.path.join(pp_dir, 'sub-*', '*', old_name)) +
                 glob.glob(os.path.join(pp_dir, 'group-*', '*', old_name)))
targets = [t for t in targets if os.path.isdir(t) and not os.path.islink(t)]
if not targets:
    sys.exit(f'No "{old_name}" folder found in {pp_dir}')

problems = []
for src in targets:
    dst = os.path.join(os.path.dirname(src), new_name)
    if os.path.exists(dst):
        problems.append(dst)
if problems:
    print('Stopped: these destination folders already exist (nothing renamed):')
    print('\n'.join('  ' + p for p in problems))
    sys.exit(1)

print(('RENAMING' if do_run else 'DRY RUN (add --run to rename)') + f': {len(targets)} folders')
for src in targets:
    dst = os.path.join(os.path.dirname(src), new_name)
    print(f'  {os.path.relpath(src, pp_dir)}  ->  {os.path.relpath(dst, pp_dir)}')
    if do_run:
        os.rename(src, dst)
print('done' if do_run else 'nothing changed')
