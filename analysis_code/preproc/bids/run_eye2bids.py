"""
-------------------------------------------------------------------------------------------------------------------------
eye2bids.py
-------------------------------------------------------------------------------------------------------------------------
Goal of the script:
Create BIDS for eyetracking folder structure. Converts edf files accordingly, copies original files to sourcedata 
folder and deletes original files. By default metedata_file.yml is expected in project_input_directory and sh script 
is saved in project_input_directory as well. 
-------------------------------------------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: project input directory path 
sys.argv[2]: project name

Option(s):
--make_copy
--delete_original
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
example: run_eye2bids.py /Users/sinakling/disks/meso_shared/  RetinoMaps  --make_copy --delete_original
-------------------------------------------------------------------------------------------------------------------------
"""

import os
import subprocess
import sys
from bids import BIDSLayout
from bids.layout import parse_file_entities
import argparse


def main(input_directory, metdadata_path, output_script_path, make_copy = True, delete_original = True):
    # Remove existing sh file 
    try:
        os.remove(f'{output_script_path}')  
    except:
        pass

    with open(output_script_path, 'w') as script_file:
        script_file.write("#!/bin/bash\n\n")
        # Make the file executable
        os.chmod(output_script_path, 0o755)
        
        # Find all .edf files
        for root, dirs, files in os.walk(input_directory):
            for file in files:
                if file.endswith(".edf"):
                    # Construct the full file path
                    edf_file_path = os.path.join(root, file)
                    # Create sourcedata folder
                    entity_dict = parse_file_entities(edf_file_path)

                    # Construct the expected output filename
                    expected_output_file = f"sub-{entity_dict['subject']}_ses-{entity_dict['session']}_task-{entity_dict['task']}_run-{entity_dict['run']}_eyetrack_recording-eye1_physio.tsv.gz"
                    
                    # Check if the output file already exists in the directory
                    if os.path.exists(os.path.join(root, expected_output_file)):
                        print(f"---Output file {expected_output_file} already exists. Skipping {edf_file_path}---")
                        continue  # Skip to the next file

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
                    # Execute conversion command
                    os.system(command)

                    if delete_original: 
                        os.remove(edf_file_path)



    print(f"Batch script generated at: {output_script_path}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert eyetracking files to BIDS format.')
    parser.add_argument('input_directory_path ', type=str, help='Project input directory path')
    parser.add_argument('input_directory_name ', type=str, help='Project input directory name')
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
    metadata_path = os.join(parent_path, 'metadata.yml')

    output_script_path = os.path.join(parent_path, 'convert_eyetracking.sh')

    # Call the main function
    main(input_directory, metadata_path, output_script_path, make_copy, delete_original)

    
 
   
