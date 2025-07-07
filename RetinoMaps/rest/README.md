# Processing of resting-state fMRI data from the RetinoMaps project using Connectome-Workbench and custom shell scripts

The workbench version used here is v2.1.0 (release date June 17, 2024). Please check the relevant installation instructions for CentOS in the unzipped folder in the meso_shared/RetinoMaps/code/ directory. 

Note you might need to update your .bashrc file to run commands from your shell with this command (and source it again):
echo 'export PATH="$PATH:/scratch/mszinte/data/RetinoMaps/code/workbench-rh_linux64-v2.1.0/bin_rh_linux64"' >> ~/.bashrc

Here we are primarily describing the steps required to:

- 1. Select the clusters from the Glasser MMP1 atlas on the fsLR 32k mesh
- 2. Downsample the task-based results from the 170k to the 91k resolution
- 3. Mask the conjunction maps from the vision pursuit and saccade tasks using the clusters by Glasser
- 4. Compute correlations using each masked task results separately for all vertices active in each cluster and averaging these values within the ROIs
- 5. Average the connectivity matrices across subjects
- 6. Prepare the files to make nice visualizations in wb_view and PyCortex:
	- Winner-take-all flatmaps
	- Pie charts

- 7. Statistical analysis (contrast maps)

Each step is described in more detail below.

## 1. Selecting the clusters from the Glasser MMP1 atlas on the fsLR 32k mesh

For this step you’ll find the scripts within the mmp1_clusters folder.

#### 1.1 Get a text file with all the label indices and rgb values

$ wb_command -cifti-label-export-table atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii INDEXMAX Glasser_labels.txt

#### 1.2 Run this bash script to match all labels you are interested in and export a txt file

$ ./filter_labels.sh

#### 1.3 Filter for the labels we are interested in

$ wb_command -cifti-label-import atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii Glasser_filtered_labels.txt atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii -discard-others;

#### 1.4 Match the filtered labels to each macro-region to generate label txt files

$ ./get_mmp1_clusters_labels_text.sh

#### 1.5 Merge labels for all the macro-regions and convert them to metric files (run the script below)

$ ./get_mmp1_clusters_metric_files.sh

#### 1.6 Run this to get the hollow seed files for visualization and stats

$ ./leave_one_out_mmp1_clusters.sh

## 2. Downsampling the task-based results from the 170k to the 91k resolution

Missing a step to describe how we obtained the left and right shape files for the filtered MMP1 ROIs.

Just run this script to perform all the steps required (see the comments in the script itself):

$ ./resample_to_91k.sh

## 3. Masking the conjunction maps using the MMP1 macro-regions

It’s also a single script that takes care of this step:

$ ./mask_task_results.sh

## 4. Computing correlations using each masked task result separately for all vertices active in each macro-region and averaging these values within the ROI bilaterally

Here we have two scripts, one for the full correlation and the other to compute the Fisher-Z of the full correlation:

$ ./compute_dtseries_corr_bilateral.sh

$ ./compute_dtseries_corr_fisher-z_bilateral.sh

## 5. Averaging the connectivity matrices across subjects

First, we need to stack all individual outputs:

$ ./cifti_stack.sh

Then, we can average across subjects using the median method:

$ ./cifti_compute_median.sh

## 6. Statistical analysis (contrast maps)

For this part, we are using the Fisher-Z outputs that are masked by the macro-regions (named hollow seed in the folders), thus excluding autocorrelation between nearby vertices. This is still a work in progress; the relevant scripts for now are:

$ ./cifti_mask_hollow_seed_stats.sh

## 7. Preparing the files to make visualizations in wb_view and PyCortex

We are trying to make at least 3 or maybe 4 kinds of plots: the individual seed functional connectivity flatmaps, the winner-take-all flatmaps, the radar plots, and the pie charts.

Some examples are below, divided by the type of plot.

Individual seed functional connectivity - Not much to say here except these are created using this script

$ ./cifti_mask_hollow_seed_viz.sh

And some examples are saved in this wb_view scene file for now: retinomaps_group.scene

Winner-take-all - These files are generated from these two scripts:

$ ./cifti_mask_hollow_seed_viz_winner-take-all.sh

After the data manipulation continues in this folder: RetinoMaps/derivatives/func_connectivity/visualizations/wb_view/wta



