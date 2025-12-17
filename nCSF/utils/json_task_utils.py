"""
-----------------------------------------------------------------------------------------
json_task_utils.py
-----------------------------------------------------------------------------------------
Goal of the script:
Utility functions for create json_task 
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

import json
import scipy.io
import os

def nCSF_jsonn(const):
    nCSF_metadata = {

        "onset": {
            "LongName": "Onset time",
            "Description": "Time at which the event occurred, relative to the start of the scan, in seconds",
            "Units": "s"
        },

        "duration": {
            "LongName": "Duration",
            "Description": "Duration of the event in seconds",
            "Units": "s"
        },

        "run_number": {
            "LongName": "Run number",
            "Description": "Indicates the number of the run"
        },

        "trial_number": {
            "LongName": "Trial number",
            "Description": "Indicates the number of the trial"
        },

        "spatial_frequency": {
            "LongName": "Spatial frequency of the stimulus",
            "Description": "Spatial frequency (cycles per degree) of the stimulus presented in the trial",
            "Units": "cycles/degree",
            "Levels": {
                "1": str(const['sf_filtCenters'][0][0]),
                "2": str(const['sf_filtCenters'][0][1]),
                "3": str(const['sf_filtCenters'][0][2]),
                "4": str(const['sf_filtCenters'][0][3]),
                "5": str(const['sf_filtCenters'][0][4]),
                "6": str(const['sf_filtCenters'][0][5]),
                "7": "break"
            }
        },

        "michelson_contrast": {
            "LongName": "Michelson contrast",
            "Description": "Contrast of the stimulus in the trial",
            "Levels": {
                "1": str(const['contValues'][0][0]),
                "2": str(const['contValues'][0][1]),
                "3": str(const['contValues'][0][2]),
                "4": str(const['contValues'][0][3]),
                "5": str(const['contValues'][0][4]),
                "6": str(const['contValues'][0][5]),
                "7": str(const['contValues'][0][6]),
                "8": str(const['contValues'][0][7]),
                "9": str(const['contValues'][0][8]),
                "10": str(const['contValues'][0][9]),
                "11": str(const['contValues'][0][10]),
                "12": str(const['contValues'][0][11]),
                "13": "break"
            }
        },

        "probe_orientation": {
            "LongName": "Orientation of the probe",
            "Description": "Orientation of the probe stimulus presented during the trial",
            "Units": "degrees",
            "Levels": {
                "1": str(const['noise_orientations'][0][0]),
                "2": str(const['noise_orientations'][0][1]),
                "3": "no probe"
            }
        },

        "response_correctness": {
            "LongName": "Correctness of response",
            "Description": "Indicates whether the participant response was correct",
            "Levels": {
                "0": "incorrect",
                "1": "correct",
                "n/a": "not applicable"
            }
        },

        "probe_time": {
            "LongName": "Probe onset time",
            "Description": "Time at which the probe appeared, relative to scan start",
            "Units": "s"
        },

        "reaction_time": {
            "LongName": "Reaction time",
            "Description": "Time from probe onset to participant response; n/a if no response",
            "Units": "s"
        }
    }

    return nCSF_metadata

def pRF_jsonn(const):
    pRF_metadata = {

        "onset": {
            "LongName": "Onset time",
            "Description": "Time at which the event occurred, relative to the start of the scan",
            "Units": "s"
        },

        "duration": {
            "LongName": "Duration",
            "Description": "Duration of the event",
            "Units": "s"
        },

        "run_number": {
            "LongName": "Run number",
            "Description": "Indicates the number of the run"
        },

        "trial_number": {
            "LongName": "Trial number",
            "Description": "Indicates the number of the trial"
        },

        "att_task": {
            "LongName": "Type of attentional task",
            "Description": "Indicates the type of attentional task performed during the run",
            "Levels": {
                "1": "task on fixation stimulus",
                "2": "task on bar stimulus"
            }
        },

        "bar_direction": {
            "LongName": "Direction of the bar movement",
            "Description": "Indicates the direction of the bar movement",
            "Units": "degrees",
            "Levels": {
                "1": "180 deg - left",
                "3": "270 deg - down",
                "5": "000 deg - right",
                "7": "090 deg - up",
                "9": "blank"
            }
        },

        "bar_period": {
            "LongName": "Bar pass period",
            "Description": "Period index of the bar pass within the run (1–9)"
        },

        "bar_step": {
            "LongName": "Bar step",
            "Description": (
                "Step number within a period: "
                "horizontal full screen (1–32), "
                "vertical full screen (1–18), "
                "aperture (1–18), or blank (1–10)"
            )
        },

        "stim_noise_ori": {
            "LongName": "Bar noise orientation",
            "Description": "Orientation of the bar noise",
            "Levels": {
                "1": "clockwise",
                "2": "counter-clockwise",
                "3": "none"
            }
        },

        "stim_stair_val": {
            "LongName": "Probe staircase kappa value",
            "Description": (
                "Kappa value of the von Mises distribution used to filter "
                "the bar orientation noise"
            ),
            "Levels": {
                "0": "noise",
                "1": str(const['noise_kappa'][0][1]),
                "2": str(const['noise_kappa'][0][2]),
                "3": str(const['noise_kappa'][0][3]),
                "4": str(const['noise_kappa'][0][4]),
                "5": str(const['noise_kappa'][0][5]),
                "6": str(const['noise_kappa'][0][6]),
                "7": str(const['noise_kappa'][0][7]),
                "8": str(const['noise_kappa'][0][8]),
                "9": str(const['noise_kappa'][0][9]),
                "10": str(const['noise_kappa'][0][10]),
                "11": str(const['noise_kappa'][0][11]),
                "12": str(const['noise_kappa'][0][12]),
                "13": str(const['noise_kappa'][0][13]),
                "14": str(const['noise_kappa'][0][14])
            }
        },

        "response_val": {
            "LongName": "Response correctness",
            "Description": "Indicates whether the participant response was correct",
            "Levels": {
                "0": "incorrect",
                "1": "correct",
                "n/a": "no response"
            }
        },

        "probe_time": {
            "LongName": "Probe onset time",
            "Description": "Time at which the oriented noise probe appeared",
            "Units": "s"
        },

        "reaction_time": {
            "LongName": "Reaction time",
            "Description": "Time from probe onset to participant response",
            "Units": "s"
        }
    }

    return pRF_metadata

def DeepMReye_jsonn(const):

    DeepMReye_metadata = {

        "onset": {
            "LongName": "Onset time",
            "Description": "Time at which the event occurred, relative to the start of the scan",
            "Units": "s"
        },

        "duration": {
            "LongName": "Duration",
            "Description": "Duration of the event",
            "Units": "s"
        },

        "run_number": {
            "LongName": "Run number",
            "Description": "Indicates the number of the run"
        },

        "trial_number": {
            "LongName": "Trial number",
            "Description": "Indicates the number of the trial"
        },

        "trial_type": {
            "LongName": "Task name",
            "Description": (
                "Task performed by the subject "
                "(intertrial, fixation, pursuit, freeview)"
            ),
            "Levels": {
                "1": str(const['task_txt'][0][0][0]),
                "2": str(const['task_txt'][0][1][0]),
                "3": str(const['task_txt'][0][2][0]),
                "4": str(const['task_txt'][0][3][0])
            }
        },

        "fixation_location": {
            "LongName": "Fixation point location",
            "Description": "Fixation location on a 5 by 5 grid spanning 18 degrees of visual angle",
            "Levels": {
                "1": str(const['fixations_postions_txt'][0][0][0]),
                "2": str(const['fixations_postions_txt'][0][1][0]),
                "3": str(const['fixations_postions_txt'][0][2][0]),
                "4": str(const['fixations_postions_txt'][0][3][0]),
                "5": str(const['fixations_postions_txt'][0][4][0]),
                "6": str(const['fixations_postions_txt'][0][5][0]),
                "7": str(const['fixations_postions_txt'][0][6][0]),
                "8": str(const['fixations_postions_txt'][0][7][0]),
                "9": str(const['fixations_postions_txt'][0][8][0]),
                "10": str(const['fixations_postions_txt'][0][9][0]),
                "11": str(const['fixations_postions_txt'][0][10][0]),
                "12": str(const['fixations_postions_txt'][0][11][0]),
                "13": str(const['fixations_postions_txt'][0][12][0]),
                "14": str(const['fixations_postions_txt'][0][13][0]),
                "15": str(const['fixations_postions_txt'][0][14][0]),
                "16": str(const['fixations_postions_txt'][0][15][0]),
                "17": str(const['fixations_postions_txt'][0][16][0]),
                "18": str(const['fixations_postions_txt'][0][17][0]),
                "19": str(const['fixations_postions_txt'][0][18][0]),
                "20": str(const['fixations_postions_txt'][0][19][0]),
                "21": str(const['fixations_postions_txt'][0][20][0]),
                "22": str(const['fixations_postions_txt'][0][21][0]),
                "23": str(const['fixations_postions_txt'][0][22][0]),
                "24": str(const['fixations_postions_txt'][0][23][0]),
                "25": str(const['fixations_postions_txt'][0][24][0])
            }
        },

        "pursuit_amplitude": {
            "LongName": "Pursuit amplitude",
            "Description": "Distance travelled by the fixation point during one trial",
            "Units": "deg",
            "Levels": {
                "1": str(const['pursuit_amps_txt'][0][0][0]),
                "2": str(const['pursuit_amps_txt'][0][1][0]),
                "3": str(const['pursuit_amps_txt'][0][2][0])
            }
        },

        "pursuit_angle": {
            "LongName": "Pursuit angle",
            "Description": "Angle of pursuit movement within an 18-degree visual field",
            "Units": "deg",
            "Levels": {
                "1": str(const['pursuit_angles_txt'][0][0][0]),
                "2": str(const['pursuit_angles_txt'][0][1][0]),
                "3": str(const['pursuit_angles_txt'][0][2][0]),
                "4": str(const['pursuit_angles_txt'][0][3][0]),
                "5": str(const['pursuit_angles_txt'][0][4][0]),
                "6": str(const['pursuit_angles_txt'][0][5][0]),
                "7": str(const['pursuit_angles_txt'][0][6][0]),
                "8": str(const['pursuit_angles_txt'][0][7][0]),
                "9": str(const['pursuit_angles_txt'][0][8][0]),
                "10": str(const['pursuit_angles_txt'][0][9][0]),
                "11": str(const['pursuit_angles_txt'][0][10][0]),
                "12": str(const['pursuit_angles_txt'][0][11][0]),
                "13": str(const['pursuit_angles_txt'][0][12][0]),
                "14": str(const['pursuit_angles_txt'][0][13][0]),
                "15": str(const['pursuit_angles_txt'][0][14][0]),
                "16": str(const['pursuit_angles_txt'][0][15][0]),
                "17": str(const['pursuit_angles_txt'][0][16][0]),
                "18": str(const['pursuit_angles_txt'][0][17][0])
            }
        },

        "image_id": {
            "LongName": "Image identifier",
            "Description": "Identifier of the image shown during the free-viewing task"
        },

    }

    return DeepMReye_metadata
    
    
def create_subject_task_events_json(base_dir, subject, session, task, run="01"):
    """
    Create the events JSON file for a given subject, session, task, and run.
    Combines screen metadata and task metadata into a single JSON file.
    """

    # Load the MATLAB file
    mat_file_fn = os.path.join(
        base_dir, "sourcedata", subject, session, "func",
        "{}_{}_task-{}_dir-PA_run-{}_matFile.mat".format(subject, session, task, run))
    mat = scipy.io.loadmat(mat_file_fn)
    config = mat["config"][0, 0]
    scr = config["scr"][0, 0]
    const = config['const'][0, 0]

    # Select the metadata function based on the task
    if task == 'pRF':
        task_metadata = pRF_jsonn(const)
    elif task == 'nCSF':
        task_metadata = nCSF_jsonn(const)
    elif task == 'DeepMReye':
        task_metadata = DeepMReye_jsonn(const)
    else:
        raise ValueError("Unknown task: {}".format(task))

    # Screen metadata
    screen_metadata = {
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
        "TaskName": task,
        "SoftwareName": "Psychtoolbox",
        "SoftwareVersion": ""
    }

    # Merge the two dictionaries into a single JSON
    full_metadata = screen_metadata.copy()  # start with screen metadata
    full_metadata.update(task_metadata)     # add task metadata

    # JSON file name
    task_json_fn = os.path.join(
        base_dir, subject, session,
        "{}_{}_task-{}_events.json".format(subject, session, task)
    )

    # Write to the JSON file
    with open(task_json_fn, "w") as f:
        json.dump(full_metadata, f, indent=4)

    print("JSON created: {}".format(task_json_fn))
        
        