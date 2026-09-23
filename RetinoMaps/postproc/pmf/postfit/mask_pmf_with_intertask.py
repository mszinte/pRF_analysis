"""
-----------------------------------------------------------------------------------------
mask_pmf_with_intertask.py
-----------------------------------------------------------------------------------------
Goal of the script:
Mask a pRF derivatives image using the vertex mask produced by make_intertask_img.py,
combining a chosen set of intertask categories
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. pmf)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python mask_pmf_with_intertask.py /scratch/mszinte/data RetinoMaps sub-01 pmf
-----------------------------------------------------------------------------------------
"""
import os
import sys
import argparse
import numpy as np
import nibabel as nb

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import make_surface_image, load_surface


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
task_names = analysis_info['analysis_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
averaging_templates = analysis_info['averaging_templates']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]

map_names_prf = analysis_info["maps_names_css"]
map_names_intertask = analysis_info['intertask_code_names']

mask_categories = ["saccade", "vision"]
mask_row_idxs = [map_names_intertask[c] for c in mask_categories]

for avg_method in avg_methods:
    for format_, extension in zip(formats, extensions):
        for task_name in task_names:
                print(f'\n{"="*72}')
                print(f'Processing: {avg_method} - {format_} - {task_name}')

                if format_ == 'fsnative':
                    # fsnative is split across two hemisphere files -- build and
                    # process each one separately
                    hemis = ['L', 'R']
                else:
                    hemis = [None]

                for hemi in hemis:
                    if format_ == 'fsnative':
                        # build paths
                        prf_fn = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/prf_derivatives/{subject}_task-{task_name}_hemi-{hemi}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_pmf-css{dm_name}_deriv.{extension}"
                        prf_img, prf_data = load_surface(fn=prf_fn)

                        intertask_fn = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/intertask/intertask_derivatives/{subject}_task-Sac-Pur-pRF_hemi-{hemi}_{preproc_prep}_{filtering}_{normalization}_loo-avg_intertask.{extension}"
                        intertask_img, intertask_data = load_surface(fn=intertask_fn)

                    else:
                        prf_fn = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/pmf/prf_derivatives/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_pmf-css{dm_name}_deriv.{extension}"
                        prf_img, prf_data = load_surface(fn=prf_fn)

                        intertask_fn = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/intertask/intertask_derivatives/{subject}_task-Sac-Pur-pRF_{preproc_prep}_{filtering}_{normalization}_loo-avg_intertask.{extension}"
                        intertask_img, intertask_data = load_surface(fn=intertask_fn)

                    if prf_data.shape[1] != intertask_data.shape[1]:
                        raise ValueError(
                            "Vertex count mismatch: pRF file has {} vertices, intertask file has {}. "
                            "Make sure both files are the same subject/hemisphere/space.".format(
                                prf_data.shape[1], intertask_data.shape[1]
                            )
                        )

                    if prf_data.shape[0] != len(map_names_prf):
                        raise ValueError(
                            "pRF file has {} maps but {} names were given ({}). "
                            "Check maps_names_css in your settings.".format(
                                prf_data.shape[0], len(map_names_prf), map_names_prf
                            )
                        )

                    # Build boolean mask
                    per_category_masks = intertask_data[mask_row_idxs, :] != 0  # (n_categories, n_vertices)

                    mask = per_category_masks.all(axis=0)

                    n_active = int(mask.sum())
                    n_total = mask.shape[0]
                    print(
                        "Mask categories {} : {}/{} vertices active ({:.1f}%)".format(
                            mask_categories, n_active, n_total,
                            100 * n_active / n_total,
                        )
                    )
                    for i, cat in enumerate(mask_categories):
                        cat_n = int(per_category_masks[i].sum())
                        print("  - '{}' alone: {}/{} vertices".format(cat, cat_n, n_total))

                    # Apply mask to every pRF map, leaving active vertices untouched (NaN)
                    masked_data = prf_data.copy().astype(float)
                    masked_data[:, ~mask] = np.nan

                    # Save
                    output_fn = prf_fn
                    os.makedirs(os.path.dirname(output_fn), exist_ok=True)
                    masked_img = make_surface_image(data=masked_data, source_img=prf_img, maps_names=map_names_prf)
                    nb.save(masked_img, output_fn)
                    print("Saved masked pRF image to {}".format(output_fn))