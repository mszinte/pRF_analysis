"""
-----------------------------------------------------------------------------------------
make_tsv_gauss.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create TSV file with all Gaussian pRF analysis output
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
sys.argv[5]: fit type ('residuals', 'bold', or 'prf')
-----------------------------------------------------------------------------------------
Output(s):
TSV file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
>> python make_tsv_gauss.py [main directory] [project name]
                            [subject num] [group] [fit_type]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python make_tsv_gauss.py /scratch/mszinte/data RetinoMaps sub-01 327 bold 
python make_tsv_gauss.py /scratch/mszinte/data RetinoMaps sub-01 327 residuals
python make_tsv_gauss.py /scratch/mszinte/data RetinoMaps sub-hcp1.6mm 327 residuals 
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import numpy as np
import pandas as pd

# Personal import
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from surface_utils import load_surface
from settings_utils import load_settings
from pycortex_utils import get_rois, set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
deriv_condition = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_names = ['SacLoc']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = ['concat']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']
maps_names_gauss_stats = analysis_info['maps_names_stats']

# Gaussian fit parameter maps
maps_names_gauss = [
    "prf_rsq",
    "prf_ecc",
    "polar_real",
    "polar_imag",
    "prf_size",
    "amplitude",
    "baseline",
    "prf_x",
    "prf_y",
    "hrf_1",
    "hrf_2",
]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define fit_type-specific filename patterns
fit_type_config = {
    'residuals': {
        'deriv_tag':  'pmf-residuals-gauss_deriv',
        'stats_tag':  'pmf-residuals-gauss_stats',
        'output_tag': 'gauss-residuals',
    },
    'bold': {
        'deriv_tag':  'pmf2-gauss_deriv',
        'stats_tag':  'pmf2-gauss_stats',
        'output_tag': 'gauss-bold',
    },
    'prf': {
        'deriv_tag':  'prf-gauss_deriv',
        'stats_tag':  'prf-gauss_stats',
        'output_tag': 'gauss-prf',
    },
}

cfg        = fit_type_config[deriv_condition]
deriv_tag  = cfg['deriv_tag']
stats_tag  = cfg['stats_tag']
output_tag = cfg['output_tag']


# Loop over averaging methods
for avg_method in avg_methods:

    # Deriv maps + stats maps combined (no pCM for gauss)
    maps_names = maps_names_gauss + maps_names_gauss_stats

    for format_, extension in zip(formats, extensions):

        # Define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

            prf_dir = "{}/{}/derivatives/pp_data/{}/{}/pmf".format(
                main_dir, project_dir, subject, format_)

            if not os.path.isdir(prf_dir):
                print(f"[SKIP] pmf_dir not found for format={format_}: {prf_dir}")
                continue

            prf_deriv_dir = "{}/pmf_derivatives".format(prf_dir)

            tsv_dir = "{}/tsv".format(prf_dir)
            os.makedirs(tsv_dir, exist_ok=True)

            for prf_task_name in prf_task_names:

                print(f'{prf_task_name} - {avg_method} - {format_}')
                tsv_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_{}.tsv'.format(
                    tsv_dir, subject, prf_task_name,
                    preproc_prep, filtering, normalization, avg_method, rois_method_format, deriv_tag)

                # Load all data
                df_rois = pd.DataFrame()

                if format_ == 'fsnative':
                    pycortex_subject = subject
                    for hemi in ['hemi-L', 'hemi-R']:

                        # Derivatives
                        deriv_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_{}.func.gii'.format(
                            prf_deriv_dir, subject, prf_task_name, hemi,
                            preproc_prep, filtering, normalization, avg_method, deriv_tag)
                        print(f'loading {deriv_fn}')
                        deriv_img, deriv_mat = load_surface(deriv_fn)

                        # Stats
                        stats_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_{}.func.gii'.format(
                            prf_deriv_dir, subject, prf_task_name, hemi,
                            preproc_prep, filtering, normalization, avg_method, stats_tag)
                        print(f'loading {stats_fn}')
                        stats_img, stats_mat = load_surface(stats_fn)

                        # Combine deriv and stats
                        all_deriv_mat = np.concatenate((deriv_mat, stats_mat))

                        # Get MMP rois img
                        roi_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(
                            main_dir, project_dir, subject, format_)
                        roi_fn = '{}/{}_{}_{}_{}_{}_rois-mmp.func.gii'.format(
                            roi_dir, subject, hemi,
                            preproc_prep, filtering, normalization)
                        roi_img, roi_mat = load_surface(roi_fn)

                        # Get MMP rois numbers tsv
                        mmp_rois_numbers_tsv_fn = os.path.join(
                            base_dir, "analysis_code", "atlas", "mmp_rois_numbers.tsv")
                        mmp_rois_numbers_df = pd.read_table(mmp_rois_numbers_tsv_fn, sep="\t")

                        # Replace roi nums by names
                        roi_mat_names = np.vectorize(
                            lambda x: dict(zip(mmp_rois_numbers_df['roi_num'],
                                               mmp_rois_numbers_df['roi_name'])).get(x, x)
                        )(roi_mat)

                        roi_verts = get_rois(subject=pycortex_subject,
                                             surf_format=format_,
                                             rois_type=rois_method_format,
                                             mask=True,
                                             rois=rois,
                                             hemis=hemi)

                        # Create and combine pandas df for each roi and hemisphere
                        print('Creating dataframe...')
                        for roi in roi_verts.keys():
                            data_dict = {
                                col: all_deriv_mat[col_idx, roi_verts[roi]]
                                for col_idx, col in enumerate(maps_names)
                            }
                            data_dict['roi'] = np.array(
                                [roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                            data_dict['roi_mmp'] = roi_mat_names[0, roi_verts[roi]]
                            data_dict['subject'] = np.array(
                                [subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                            data_dict['hemi'] = np.array(
                                [hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                            data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                            df_rois = pd.concat(
                                [df_rois, pd.DataFrame(data_dict)], ignore_index=True)

                elif format_ == '170k':
                    pycortex_subject = pycortex_subject_template

                    # Derivatives
                    deriv_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}.dtseries.nii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, deriv_tag)
                    print(f'loading {deriv_fn}')
                    deriv_img, deriv_mat = load_surface(deriv_fn)

                    # Stats
                    stats_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}.dtseries.nii'.format(
                        prf_deriv_dir, subject, prf_task_name,
                        preproc_prep, filtering, normalization, avg_method, stats_tag)
                    print(f'loading {stats_fn}')
                    stats_img, stats_mat = load_surface(stats_fn)

                    # Combine deriv and stats
                    all_deriv_mat = np.concatenate((deriv_mat, stats_mat))

                    # Get MMP rois img
                    roi_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(
                        main_dir, project_dir, subject, format_)
                    roi_fn = '{}/{}_{}_{}_{}_rois-mmp.dtseries.nii'.format(
                        roi_dir, subject,
                        preproc_prep, filtering, normalization)
                    roi_img, roi_mat = load_surface(roi_fn)

                    # Get MMP rois numbers tsv
                    mmp_rois_numbers_tsv_fn = os.path.join(
                        base_dir, "analysis_code", "atlas", "mmp_rois_numbers.tsv")
                    mmp_rois_numbers_df = pd.read_table(mmp_rois_numbers_tsv_fn, sep="\t")

                    # Replace roi nums by names
                    roi_mat_names = np.vectorize(
                        lambda x: dict(zip(mmp_rois_numbers_df['roi_num'],
                                           mmp_rois_numbers_df['roi_name'])).get(x, x)
                    )(roi_mat)

                    # Create and combine pandas df for each roi and hemisphere
                    print('Creating dataframe...')
                    for hemi in ['hemi-L', 'hemi-R']:
                        roi_verts = get_rois(subject=pycortex_subject,
                                             surf_format=format_,
                                             rois_type=rois_method_format,
                                             mask=True,
                                             rois=rois,
                                             hemis=hemi)

                        for roi in roi_verts.keys():
                            data_dict = {
                                col: all_deriv_mat[col_idx, roi_verts[roi]]
                                for col_idx, col in enumerate(maps_names)
                            }
                            data_dict['roi'] = np.array(
                                [roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                            data_dict['roi_mmp'] = roi_mat_names[0, roi_verts[roi]]
                            data_dict['subject'] = np.array(
                                [subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                            data_dict['hemi'] = np.array(
                                [hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                            data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                            df_rois = pd.concat(
                                [df_rois, pd.DataFrame(data_dict)], ignore_index=True)

                # Fill NaN for HRF columns (may not be fitted in all cases)
                df_rois[['hrf_1', 'hrf_2']] = df_rois[['hrf_1', 'hrf_2']].fillna('non_computed')
                print('Saving tsv: {}'.format(tsv_fn))
                df_rois.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))