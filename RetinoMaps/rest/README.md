# Processing of resting-state fMRI data from the RetinoMaps project

Most of the processing is carried out using Connectome-Workbench. The version used here is v2.1.0 (release date June 17, 2024). Please check the relevant installation instructions for CentOS in the unzipped folder in the meso_shared/RetinoMaps/code/ directory. 

Note you might need to update your .bashrc file to run commands from your shell with this command (and source it again):
echo 'export PATH="$PATH:/scratch/mszinte/data/RetinoMaps/code/workbench-rh_linux64-v2.1.0/bin_rh_linux64"' >> ~/.bashrc
Every time you log into the meso source the file (source .bashrc), if any script doesn't run correctly.

The rest of the pipeline uses bash scripts and Nilearn to compute the partial correlations on the surface.

The primary steps of the pipeline are:

- 1. Selecting the clusters from the MMP1 atlas (Glasser et al. 2016) on the fsLR 32k mesh
- 2. Downsampling the task-based results from the 170k to the 91k fsLR resolution
- 3. Masking the conjunction maps from the vision pursuit and saccade tasks using the MMP1 clusters
- 4. Computing full and partial correlations using each masked task result
- 5. Averaging the connectivity matrices across subjects
- 6. Preparing the files to make nice visualizations in wb_view and PyCortex:
	- Winner-take-all flatmaps
	- Violin plots of cluster's averaged r values

Each step is described in more detail below.

### 1. Selecting the clusters from the Glasser MMP1 atlas on the fsLR 32k mesh

For this step, you’ll find the scripts within the mmp1_clusters folder.

##### 1.1 Get a text file with all the label indices and RGB values
```bash
$ wb_command -cifti-label-export-table atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii INDEXMAX Glasser_labels.txt
```
#### 1.2 Run this bash script to match all labels you are interested in and export a txt file
```bash
$ ./filter_labels.sh
```
#### 1.3 Filter for the labels we are interested in
```bash
$ wb_command -cifti-label-import atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii Glasser_filtered_labels.txt atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii -discard-others;
```
#### 1.4 Match the filtered labels to each macro-region to generate label txt files
```bash
$ ./get_mmp1_clusters_labels_text.sh
```
#### 1.5 Merge labels for all the macro-regions and convert them to metric files (run the script below)
```bash
$ ./get_mmp1_clusters_metric_files.sh
```
#### 1.6 Run this to get the hollow seed files for visualization and stats
```bash
$ ./leave_one_out_mmp1_clusters.sh
```
### 2. Downsampling the task-based results from the 170k to the 91k resolution of the fsLR template

Missing a step to describe how we obtained the left and right shape files for the filtered MMP1 ROIs.

Just run this script to perform all the steps required (see the comments in the script itself):
```bash
$ ./resample_to_91k.sh
```
### 3. Masking the conjunction maps using the MMP1 macro-regions (by cluster or parcel)

Check the two scripts in the masking folder:
```bash
$ ./mask_task_results_by_cluster.sh

$ ./mask_task_results_by_parcel.sh
```
### 4. Computing correlations using each masked task result separately for all vertices active in each cluster and averaging these values bilaterally

Here we have two scripts, one for the full correlation and the other to compute the Fisher-Z transformed values of the full correlation:
```bash
$ ./compute_dtseries_corr_bilateral.sh

$ ./compute_dtseries_corr_fisher-z_bilateral.sh
```
Computing the partial correlations using Nilearn (importantly, the targets are defined by the MMP parcellation and not the task-defined vertices):

nilearn_partial_corr_cluster-task_by_mmp-parcel.py

### 5. Averaging the connectivity matrices across subjects

First, we need to stack all individual outputs:
```bash
$ ./cifti_stack.sh
```
Then, we can average across subjects using the median method:
```bash
$ ./cifti_compute_median.sh
```
### 6. Preparing the files to make visualizations in wb_view and eventually, PyCortex

We focus on two main plots: the winner-take-all flatmaps for full and partial correlation parcellated results, and the violin plots with averaged r values within the MMP1 clusters.

Some examples are below, divided by the type of plot.

Winner-take-all - These files are generated from these two scripts:
```bash
$ ./cifti_mask_hollow_seed_viz_winner-take-all.sh
```

Some examples are saved in this wb_view scene file for now: retinomaps_group.scene

After the data manipulation continues in this folder: meso_shared/RetinoMaps/derivatives/pp_data/group/91k/visualizations



