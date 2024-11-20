"""
-------------------------------------------------------------------------------------------------------------------------
run_eye2bids.py
-------------------------------------------------------------------------------------------------------------------------
Goal of the script:
Create BIDS for eyetracking folder structure. Converts edf files accordingly, copies original files to sourcedata 
folder and deletes original files. By default metedata_file.yml is expected in project_input_directory and sh script 
is saved in project_input_directory as well. 
-------------------------------------------------------------------------------------------------------------------------
Input(s):
project input directory path 
project name

Option(s):
--make_copy : copy the files from the BIDS folder to the sourcedata folder
--delete_original: delete original files from the BIDS folder
-------------------------------------------------------------------------------------------------------------------------
Output(s):
BIDS converted eyetracking file structure: 
sub-02_ses-01_task-pRF_run-03_eyeData.edf -> sub-02_ses-01_task-pRF_run-03_eyeData_recording-eye1_physio.json 
                                             sub-02_ses-01_task-pRF_run-03_eyeData_recording-eye1_physio.tsv.gz
                                             sub-02_ses-01_task-pRF_run-03_eyeData_recording-eye1_physioevents.json
                                             sub-02_ses-01_task-pRF_run-03_eyeData_recording-eye1_physiosevents.tsv.gz
-------------------------------------------------------------------------------------------------------------------------
To run locally: 
cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
python run_eye2bids.py /path/to/project project_name --make_copy --delete_original  
example: python run_eye2bids.py /Users/sinakling/disks/meso_shared/  RetinoMaps  --make_copy --delete_original
-------------------------------------------------------------------------------------------------------------------------
"""


# General imports
import os
import subprocess
import sys
import ipdb
import re
deb=ipdb.set_trace

# Import BIDS if available
try:
    from bids import BIDSLayout
    from bids.layout import parse_file_entities
except ImportError:
    print("Error: The 'pybids' package is required. Install it with 'pip install pybids'.")
    sys.exit(1)
import argparse


def main(input_directory, metadata_path, output_script_path, make_copy=True, delete_original=True):
    # Check if input directory exits 
    if not os.path.exists(input_directory):
        print(f"Error: Input directory '{input_directory}' does not exist. Please provide a valid path.")
        sys.exit(1)
    # Check if the metadata file exists
    if not os.path.exists(metadata_path):
        print(f"Error: Metadata file not found at {metadata_path}. Exiting script.")
        sys.exit(1)  # Exit the script if metadata file is missing

    # Remove existing sh file 
    try:
        os.remove(f'{output_script_path}')  
    except:
        pass

    with open(output_script_path, 'w') as script_file:
        script_file.write("#!/bin/bash\n\n")
        # Make the file executable
        os.chmod(output_script_path, 0o755)
        
        print("Looking for edf files...")
        
        # Find all .edf files
        for root, dirs, files in os.walk(input_directory):

            # Create the new filename by removing '_eyedata'
            for file in files:
                if '_eyeData' in file:
                    new_file_name = file.replace('_eyeData', '')
                    old_file_path = os.path.join(root, file)
                    new_file_path = os.path.join(root, new_file_name)
                    print(f'Renaming {old_file_path} to {new_file_path}')
                    os.rename(old_file_path, new_file_path)
            
            if 'sourcedata' in dirs:
                dirs.remove('sourcedata')
            
            for file in files:

                # Check if any output file matching the pattern already exists in the directory
                eye_pattern = re.compile(r"sub-.*_ses-.*_task-.*_run-.*_recording-.*_physio\.tsv\.gz")
                if eye_pattern.match(file):
                    print(f"---Output file {file} already exists. Skipping {file}---")
                    continue  # Skip to the next file
                
                if file.endswith(".edf"):
                    # Construct the full file path
                    edf_file_path = os.path.join(root, file)
                    
                    # Create sourcedata folder
                    entity_dict = parse_file_entities(edf_file_path)

                    # Create Sourcedata directory
                    if make_copy: 
                        sourcedata_folder = '{}/sourcedata/sub-{}/ses-{}/{}'.format(input_directory, 
                                                                            entity_dict['subject'], 
                                                                            entity_dict['session'], 
                                                                            entity_dict['datatype'])

                        os.makedirs('{}'.format(sourcedata_folder), exist_ok=True)

                        # Copy edf files into created sourcedata folder 
                        os.system(f'cp {edf_file_path} {sourcedata_folder}')
                        print(f"---Succesfully copied {edf_file_path}---")

                    # Construct the command and write it to the .sh file 
                    command = f"eye2bids --input_file {edf_file_path} --metadata_file {metadata_path} --output_dir {root}\n"  
                    script_file.write(command)
                    os.system(command)

                    # Delete original files and edf2bids .asc files
                    if delete_original:
                        print("Removing original files and .asc files")
                        os.remove(edf_file_path)
                        base_name = os.path.splitext(edf_file_path)[0]
                        events_file_path = f"{base_name}_events.asc"
                        samples_file_path = f"{base_name}_samples.asc"
                        for file_path in [events_file_path, samples_file_path]:
                            if os.path.exists(file_path):
                                os.remove(file_path)

    print(f"Batch script generated at: {output_script_path}")

# Runner
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert eyetracking files to BIDS format.')
    parser.add_argument('input_directory_path', type=str, help='Project input directory path')
    parser.add_argument('input_directory_name', type=str, help='Project input directory name')
    parser.add_argument('--make_copy', action='store_true', help='Make a copy of the original files')
    parser.add_argument('--delete_original', action='store_true', help='Delete original files after copying')
    
    args = parser.parse_args()

    # Get inputs
    input_directory_path = args.input_directory_path
    input_directory_name = args.input_directory_name
    input_directory = os.path.join(input_directory_path, input_directory_name)
    make_copy = args.make_copy
    delete_original = args.delete_original
    parent_path = os.path.abspath(os.path.join(os.getcwd(), "../../.."))
    name_path = os.path.join(parent_path,  input_directory_name)  
    metadata_path = os.path.join(name_path, 'metadata.yml')   # Required to have metadata.yml file in path: /../../projects/pRF_analysis/project_name/metadata.yml
    output_script_path = os.path.join(name_path, 'eyetracking/convert_eyetracking.sh')

    # Call the main function
    main(input_directory, metadata_path, output_script_path, make_copy, delete_original)
