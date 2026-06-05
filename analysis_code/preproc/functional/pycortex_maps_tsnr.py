"""
-----------------------------------------------------------------------------------------
pycortex_maps_tsnr.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots of tSNR maps
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save map in svg (y/n)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/functional/
2. run python command
>> python pycortex_maps_tsnr.py [main directory] [project name] [subject num] [save_svg]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/functional/
python pycortex_maps_tsnr.py ~/disks/meso_S/data amblyo7T_prf sub-01 n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug import
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import glob
import cortex
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from settings_utils import load_settings
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg = sys.argv[4]

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
tsnr_method = analysis_info['tsnr_method']
tsnr_scale = analysis_info['flatmap_tsnr_scale']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Colormap
cmap_tsnr = 'viridis'

# Loop over formats
for format_, extension, pycortex_subject in zip(
        formats, extensions, [subject, pycortex_subject_template]):

    for method in (['standard', 'robust'] if tsnr_method == 'both' else [tsnr_method]):
        fn_label = 'tSNR' if method == 'standard' else 'tSNR-robust'

        tsnr_dir = "{}/{}/derivatives/pp_data/{}/{}/func/{}_{}".format(
            main_dir, project_dir, subject, format_, preproc_prep, filtering)

        flatmaps_dir = '{}/{}/derivatives/pp_data/{}/{}/tsnr/pycortex/flatmaps_tsnr'.format(
            main_dir, project_dir, subject, format_)
        datasets_dir = '{}/{}/derivatives/pp_data/{}/{}/tsnr/pycortex/datasets_tsnr'.format(
            main_dir, project_dir, subject, format_)
        os.makedirs(flatmaps_dir, exist_ok=True)
        os.makedirs(datasets_dir, exist_ok=True)

        volumes = {}

        if format_ == 'fsnative':
            tsnr_fns_L = sorted(glob.glob(
                "{}/{}_*hemi-L*{}.func.gii".format(tsnr_dir, subject, fn_label)))
            tsnr_fns_L += sorted(glob.glob(
                "{}/{}_*{}_avg.func.gii".format(tsnr_dir, subject, fn_label)))
            
            if not tsnr_fns_L:
                print('[SKIP] No {} files found for format={}'.format(fn_label, format_))
                continue

            
            for fn_L in tsnr_fns_L:
                fn_R = fn_L.replace('hemi-L', 'hemi-R')
                if not os.path.exists(fn_R):
                    print('[SKIP] Missing hemi-R for {}'.format(fn_L))
                    continue

                map_name = os.path.basename(fn_L).replace('.func.gii', '').replace(
                    '{}_'.format(subject), '').replace('hemi-L_', '')
                
                print('Processing: {}'.format(map_name))

                results = load_surface_pycortex(L_fn=fn_L, R_fn=fn_R)
                tsnr_data = results['data_concat'][0, :]

                param_tsnr = {'data': tsnr_data,
                              'alpha': tsnr_data * 0 + 1,
                              'cmap': cmap_tsnr,
                              'vmin': tsnr_scale[0],
                              'vmax': tsnr_scale[1],
                              'cbar': 'discrete',
                              'cortex_type': 'VertexRGB',
                              'description': '{}: {}'.format(fn_label, map_name),
                              'curv_brightness': 0.1,
                              'curv_contrast': 0.25,
                              'add_roi': save_svg,
                              'cbar_label': 'tSNR',
                              'with_labels': True,
                              'with_borders': True,
                              'overlay_fn': 'overlays_rois-mmp.svg',
                              'subject': pycortex_subject,
                              'roi_name': map_name}

                volume_tsnr = draw_cortex(**param_tsnr)
                plt.savefig('{}/{}_{}.pdf'.format(flatmaps_dir, subject, map_name))
                plt.close()
                volumes[param_tsnr['description']] = volume_tsnr

        elif format_ == '170k':
            tsnr_fns = sorted(glob.glob(
                "{}/{}_*{}.dtseries.nii".format(tsnr_dir, subject, fn_label)))
            tsnr_fns += sorted(glob.glob(
                "{}/{}_*{}_avg.dtseries.nii".format(tsnr_dir, subject, fn_label)))
            if not tsnr_fns:
                print('[SKIP] No {} files found for format={}'.format(fn_label, format_))
                continue

            for fn in tsnr_fns:
                map_name = os.path.basename(fn).replace('.dtseries.nii', '').replace(
                    '{}_'.format(subject), '')
                print('Processing: {}'.format(map_name))

                results = load_surface_pycortex(brain_fn=fn)
                tsnr_data = results['data_concat'][0, :]

                param_tsnr = {'data': tsnr_data,
                              'alpha': tsnr_data * 0 + 1,
                              'cmap': cmap_tsnr,
                              'vmin': tsnr_scale[0],
                              'vmax': tsnr_scale[1],
                              'cbar': 'discrete',
                              'cortex_type': 'VertexRGB',
                              'description': '{}: {}'.format(fn_label, map_name),
                              'curv_brightness': 0.1,
                              'curv_contrast': 0.25,
                              'add_roi': save_svg,
                              'cbar_label': 'tSNR',
                              'with_labels': True,
                              'with_borders': True,
                              'overlay_fn': 'overlays_rois-mmp.svg',
                              'subject': pycortex_subject,
                              'roi_name': map_name}

                volume_tsnr = draw_cortex(**param_tsnr)
                plt.savefig('{}/{}_{}.pdf'.format(flatmaps_dir, subject, map_name))
                plt.close()
                volumes[param_tsnr['description']] = volume_tsnr

        # Save single HDF for all maps of this format/method
        if volumes:
            dataset_file = "{}/{}_{}.hdf".format(datasets_dir, subject, fn_label)
            if os.path.exists(dataset_file):
                os.system("rm -fv {}".format(dataset_file))
            dataset = cortex.Dataset(data=volumes)
            dataset.save(dataset_file)
            print('Saved dataset: {}'.format(dataset_file))