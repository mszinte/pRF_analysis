#!/usr/bin/env python3
import os
from pathlib import Path

main_folder = '/scratch/mszinte/data/amblyo7T_prf/derivatives/fmriprep/fmriprep'
old_prefix = 'sub-03'
new_prefix = 'sub-10'

# Find all matching files
files_to_rename = []
for root, dirs, files in os.walk(main_folder):
    # Skip derivatives folder
    dirs[:] = [d for d in dirs if d != 'derivatives']
    
    # Limit depth to 2 levels
    depth = root.replace(main_folder, '').count(os.sep)
    if depth > 2:
        dirs[:] = []
        continue
    
    for file in files:
        if file.startswith(old_prefix):
            old_path = Path(root) / file
            new_name = file.replace(old_prefix, new_prefix, 1)
            new_path = Path(root) / new_name
            files_to_rename.append((old_path, new_path))

# Show preview
print("\nFiles to be renamed:")
print("=" * 50)
for old, new in files_to_rename:
    print(f"{old} -> {new}")

# Confirm
print("\n" + "=" * 50)
answer = input("Proceed with renaming? (y/n): ")

if answer.lower() == 'y':
    for old, new in files_to_rename:
        old.rename(new)
    print("Renaming complete!")
else:
    print("Cancelled.")