# Processing of resting-state fMRI data from the RetinoMaps project

Most of the processing is carried out using Connectome-Workbench. The version used here is v2.1.0 (release date June 17, 2024). Please check the relevant installation instructions for CentOS in the unzipped folder in the meso_shared/RetinoMaps/code/ directory. 

Note you might need to update your .bashrc file to run commands from your shell with this command (and source it again):
echo 'export PATH="$PATH:/scratch/mszinte/data/RetinoMaps/code/workbench-rh_linux64-v2.1.0/bin_rh_linux64"' >> ~/.bashrc
Every time you log into the meso source the file (source .bashrc), if any script doesn't run correctly.

The rest of the pipeline uses bash scripts and Nilearn to compute the partial correlations on the surface.
Please don't forget to activate your environment with conda activate pRF_env.

The primary steps of the pipeline are:

- 1. Selecting the macro-regions from the MMP1 atlas (Glasser et al. 2016) on the fsLR 32k mesh
- 2. Downsampling the task-based results from the 170k to the 91k fsLR resolution
- 3. Masking the conjunction maps from the vision, saccade and pursuit tasks using the MMP1 clusters
- 4. Computing full and partial correlations using each masked task result intra-hemispherically
- 5. Averaging the 12 seed x 53 targets connectivity matrices across subjects
- 6. Preparing the files to make nice visualizations in wb_view and PyCortex:
	- Winner-take-all flat maps by subject and at the group-level
	- Violin plots of macro-regions's averaged r values by hemisphere (stacked)

Each step is described in more detail below.

### 1. Selecting the macro-regions from the Glasser MMP1 atlas on the fsLR 32k mesh

For this step, you’ll find the scripts within the atlas folder.

##### 1.1 Get a text file with all the label indices and RGB values
```bash
$ wb_command -cifti-label-export-table /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii INDEXMAX /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1/mmp1_labels.txt
```
#### 1.2 Go to the atlas folder and run this bash script to match all labels you are interested in and export a txt file
```bash
$ ./filter_labels.sh
```
#### 1.3 Filter for the labels we are interested in
```bash
$ wb_command -cifti-label-import /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1/Glasser_filtered_labels.txt /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii -discard-others;
```
#### 1.4 Match the filtered labels to each macro-region to generate label txt files
```bash
$ ./get_mmp1_macro-regions_labels_text.sh
```
#### 1.5 Merge labels for all the macro-regions and convert them to metric files
```bash
$ ./get_mmp1_macro-regions_metric_files.sh
```
#### 1.6 Run this to get the hollow seed files for visualization and stats
```bash
$ ./leave_one_out_mmp1_macro-regions.sh
```
### 2. Downsampling the tasks results from the 170k to the 91k resolution of the fsLR template

Run this script to perform all the steps required (see the comments in the script itself):
```bash
$ ./resample_to_fsLR91k_adap-bary.sh
```
### 3. Masking the conjunction maps using the MMP1 macro-regions

Check the script in the masking folder:
```bash
$ ./mask_task_results_by_macro-region.sh
```
### 4. Computing functional connectivity using each macro-region separately by hemisphere (and per run)

The relevant scripts live in the correlation folder.

Full correlation are computed with connectome-workbench to obtain outputs as Pearson r and Fisher-z transformed values:
```bash
$ ./connectome-workbench_seed-task_by_mmp-parcel_full_corr_by_hemi.sh
```
Partial correlations are computed using Nilearn mirroring the workbench approach (except timeseries are averaged within macro-regions seeds and targets):
```bash
$ python nilearn_partial_corr_seed-task_by_mmp-parcel_by_hemi.py
```

### 5. Averaging the parcellated connectivity matrices across subjects:

Two scripts in the stats folder that by default run on fisher-z outputs:
```bash
$ python group_stats_full_corr.py
$ python group_stats_partial_corr.py
```

### 6. Preparing the files to make visualizations in PyCortex

In the visualization folder, two scripts handle the winner-take-all subject-wise and per run computation for full and partial correlation parcellated results.
```bash
$ python wta_full_corr_by_subject_by_hemi_per_run.py
$ python wta_partial_corr_by_subject_by_hemi_per_run.py
```

For inspecting these results using wb_view, you can run the following bash scripts to produce the appropriate dlabel files:
```bash
$ ./generate_full_corr_wta_dlabel_file.sh
$ ./generate_partial_corr_wta_dlabel_file.sh
```
Some example visualizations are going to be released soon as a .scene file for workbench view.

Violin plots are instead created using these two scripts:
```bash
$ python violin_plots_workbench_full_corr_seeds_by_target.py
$ python violin_plots_nilearn_partial_corr_seeds_by_target.py
```