# About
---
*We here study cortical areas involved in both vision and eye movements processes in 20 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** and **HCP 170k or 91k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

# Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

# Data stucture and BIDS
---
- [x] Copy relevant data from PredictEye [copy_data.py](preproc/bids_copy_data.sh) 
- [x] Change the 'task' to 'task_condition' column name in event.tsv files to avoid BIDS problems [correct_events_files.ipynb](preproc/correct_events_files.ipynb)
- [x] Deface participants t1w image [bidsonym_sbatch.py](../analysis_code/preproc/bids/bidsonym_sbatch.py) 
- [x] BIDS eyetracking conversion [run_eye2bids.py](../analysis_code/preproc/bids/run_eye2bids.py)
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]

# Structural preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Create sagittal view video before manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Manual edit of brain segmentation [pial_edits.sh](../analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](../analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] Create sagittal view video after manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Flattening of the cortex using [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py) (with or without manual cut [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh))

# Functional preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, anat [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py) or ro run on server [preproc_end_sbatch.py](../analysis_code/preproc/functional/preproc_end_sbatch.py)
- [x] Averaging across runs [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](preproc/functional/averaging.py)

# Inter-run correlations 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] Compute inter-run correlation [compute_run_corr_sbatch](../analysis_code/preproc/functional/compute_run_corr_sbatch.py)
- [x] Make maps with pycortex [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](../analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**170k**) in which individual results are averaged.</br> 

- [x] Compute inter-run correlation for **template_avg** [compute_run_corr.py](../analysis_code/preproc/functional/compute_run_corr.py)
- [x] Make maps with pycortex for **hcp1.6mm** [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py)

# pRF 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

#### PRF Gaussian fit
- [x] Create the visual matrix design [vdm_builder_prf.py](postproc/prf/fit/vdm_builder_prf.py)
- [x] Run pRF gaussian fit [prf_submit_gaussfit_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_gauss_jobs.py)
- [x] Compute pRF gaussian fit derivatives [compute_gauss_derivatives.py](../analysis_code/postproc/prf/postfit/compute_gauss_derivatives.py)
- [x] Make pRF maps with pycortex [pycortex_maps_gauss.py](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.py) or [pycortex_maps_gauss.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.sh)

#### PRF ROIs
- [x] Create 170k MMP rois masks [create_rois-mmp_npz.py](../analysis_code/atlas/create_rois-mmp_npz.py)
- [x] Draw individual ROI on fsnative data using Inkscape and create masks rois mask [create_rois-drawn_npz.py](../analysis_code/postproc/prf/postfit/create_rois-drawn_npz.py)
- [x] Project MMP atlas on fsnative surface using freesurfer [freesurfer_project_mmp_fsnative.py](../analysis_code/postproc/prf/postfit/freesurfer_project_mmp_fsnative.py)
- [x] Make fsnative roi mmp npz and overlays [fsnative_mmp_rois.py](../analysis_code/postproc/prf/postfit/fsnative_mmp_rois.py)
- [x] Make ROIS files [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

#### PRF CSS fit
- [x] CSS fit within the ROIs [prf_submit_css_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [css_stats_sbatch.py](../analysis_code/postproc/prf/postfit/css_stats_sbatch.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](../analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make maps with pycortex [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [x] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) or [make_rois_fig_tsv.sh](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.sh)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](../analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [x] Merge all figures [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**170k**) in which individual results are averaged and on an ROI-based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

#### PRF Gaussian fit
- [x] Compute pRF gaussian grid fit derivatives for **template_avg** [compute_gauss_derivatives.py](../analysis_code/postproc/prf/postfit/compute_gauss_derivatives.py)
- [x] Make pRF maps with pycortex for **hcp1.6mm**  [pycortex_maps_gauss.py](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.py)

#### PRF ROIs
- [x] Make ROIS files for **template_avg** [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs for **hcp1.6mm** [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py)

#### PRF CSS fit
- [x] Compute CSS statistics for **template_avg** [compute_css_stats.py](../analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives for **template_avg** [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) for **template_avg** [compute_css_pcm.py](../analysis_code/postproc/prf/postfit/compute_css_pcm.py)
- [x] Make maps with pycortex for **hcp1.6mm** [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics and CM for **hcp1.6mm** [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM for **hcp1.6mm** and **group** [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM for **hcp1.6mm** and **group** [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [x] Merge all figures for **hcp1.6mm** and **group** [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)

# GLM 
--- 
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] Run GLM for the differents tasks [glm_fit.py](glm/fit/glm_fit.py) or [glm_sbatch.py](glm/fit/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](glm/postfit/compute_glm_stats.py) or [glm_stats_sbatch.py](glm/postfit/glm_stats_sbatch.py)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**sub-170k**) in which individual results are averaged and on an ROI-based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

- [x] Compute GLM statistics for **sub-170k** [compute_glm_stats.py](glm/postfit/compute_glm_stats.py) or [glm_stats_sbatch.py](glm/postfit/glm_stats_sbatch.py)

# Inter task analysis 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] Make intertasks image [make_intertask_img.py](intertask/make_intertask_img.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_tsv.py](intertask/make_intertask_tsv.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics, CM  and GLM results [make_intertask_rois_fig_tsv_sbatch.py](intertask/make_intertask_rois_fig_tsv_sbatch.py) 
- [x] Make figure specific TSV with GLM results for active vertex [make_active_vert_fig_tsv.py](intertask/make_active_vert_fig_tsv.py) 
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_rois_fig.py](intertask/make_intertask_rois_fig.py) or [make_intertask_rois_fig.sh](intertask/make_intertask_rois_fig.sh)
- [x] Make final statistical maps maps with pycortex [pycortex_maps_intertask.py](intertask/pycortex_maps_intertask.py) or [pycortex_maps_intertask.sh](intertask/pycortex_maps_intertask.sh)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**sub-170k**) in which individual results are averaged and on an ROI-based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

- [x] Make intertasks image for **sub-170k** [make_intertask_img.py](intertask/make_intertask_img.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics, CM and GLM results for **sub-170k** [make_intertask_tsv.py](intertask/make_intertask_tsv.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics, CM  and GLM results for **sub-170k** and **group** [make_intertask_rois_fig_tsv.py](intertask/make_intertask_rois_fig_tsv.py) or [make_intertask_rois_fig_tsv_sbatch.py](intertask/make_intertask_rois_fig_tsv_sbatch.py)
- [x] Make figure specific TSV with GLM results for active vertex for **sub-170k** and **group** [make_active_vert_fig_tsv.py](intertask/make_active_vert_fig_tsv.py) 
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics, CM and GLM results for **sub-170k** and **group** [make_intertask_rois_fig.py](intertask/make_intertask_rois_fig.py) or [make_intertask_rois_fig.sh](intertask/make_intertask_rois_fig.sh)
- [x] Make intertask statistical maps maps with pycortex  for **sub-170k** [pycortex_maps_intertask.py](intertask/pycortex_maps_intertask.py) or [pycortex_maps_intertask.sh](intertask/pycortex_maps_intertask.sh)

# Resting-state preprocessing
---
- [x] Preprocess the resting-state data and output it in the fsLR-91k resolution [fmriprep_sbatch_ica-aroma.py](https://github.com/mszinte/pRF_analysis/blob/main/analysis_code/preproc/functional/fmriprep_sbatch_ica-aroma.py)

# Resting-state postprocessing
---
- [x] Extract motion components with ICA-AROMA using fMRIPost-AROMA [fmripost_aroma_sbatch.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/rest/fmripost_aroma_sbatch.py)
- [x] Post-process and denoise the data using XCP-D [xcp-d_aroma_sbatch.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/rest/xcp-d_aroma_sbatch.py)

# Resting-state analysis
---
- [x] Compute seed-based (from the tasks' conjunction results) functional connectivity on the dense timeseries with connectome workbench: Pearson correlation ([compute_dtseries_corr_bilateral.sh](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/rest/correlations/) and Fisher-z transformed Pearson correlation ([compute_dtseries_corr_fisher-z_bilateral.sh](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/rest/correlations/compute_dtseries_corr_fisher-z_bilateral.sh)
- [ ] Compute seed-based partial correlation
- [ ] Compute winner-take-all results

# Resting-state results visualization
---
TO DO

# WEBGL
---
https://invibe.nohost.me/predicteye/

- [x] Make subject WEBGL with pycortex [pycortex_webgl_css.py](webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](../analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](webgl/publish_webgl.py)

# Eyetracking preprocessing 
---
- [x] Generate experimental design matrix [create_eye_tracking_DM.py](../analysis_code/preproc/bids/create_eye_tracking_DM.py)
- [ ] Eyetrack preprocessing [eyetrack_preproc.py](eyetracking/eyetrack_preproc.py)

# Eyetracking postprocessing
---
- [ ] Extract trigger timestamps [extract_triggers.py](eyetracking/extract_triggers.py)
- [ ] Extract saccades [extract_saccades.py](eyetracking/extract_saccades.py)
- [ ] Create saccade model prediction timeseries [saccade_prediction.py](eyetracking/saccade_prediction.py)
- [ ] Generate predictions timeseries and extract stats [generate_prediction.py](eyetracking/generate_prediction.py)
- [ ] Create timeseries figures [timeseries_figures.py](eyetracking/timeseries_figures.py)
- [ ] Create stats figure [stats_figures.py](eyetracking/stats_figures.py)

# pRF behaviour
---
- [x] Make pRF behaviour figure TSV [make_prf_beh_fig_tsv.py](pRF_beh/make_prf_beh_fig_tsv.py)
- [x] Make pRF behaviour figure [make_prf_beh_fig.py](pRF_beh/make_prf_beh_fig.py)


# pMF 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

#### PMF Gaussian fit
- [x] Create the visual matrix design with eye movements (retinal view vdm) [vdm_builder_sacloc.py](postproc/pmf/fit/vdm_builder_sacloc.py)
- [ ] Concatonate SacLoc runs [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](preproc/functional/averaging.py)