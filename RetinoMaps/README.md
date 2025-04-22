# RetinoMaps
## About
---
*We here study cortical areas involved in both vision and eye movements in 20 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** and **HCP 170k or 91k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

---
## Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

## To Do 
---
- [ ] take care of glm analysis
- [ ] add quality check file

## Data analysis
---
### Preprocessing
- [x] Copy relevant data from PredictEye [copy_data.py](preproc/bids_copy_data.sh) 
- [x] Change the 'task' to 'task_condition' column name in event.tsv files to avoid BIDS problems [correct_events_files.ipynb](preproc/correct_events_files.ipynb)

### pRF behaviour analysis
- [x] Make pRF behaviour figure TSV [make_prf_beh_fig_tsv.py](pRF_beh/make_prf_beh_fig_tsv.py)
- [x] Make pRF behaviour figure [make_prf_beh_fig.py](pRF_beh/make_prf_beh_fig.py)

### Make project WEBGL
- [x] Make subject WEBGL with pycortex [pycortex_webgl_css.py](webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](webgl/publish_webgl.py)

### GLM analysis
- [x] Run GLM for the differents tasks [glm_sbatch.py](glm/fit/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](glm/postfit/compute_glm_stats.py)

### Inter task analysis
- [x] Make intertasks image [make_intertask_img.py](intertask/make_intertask_img.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_tsv.py](intertask/make_intertask_tsv.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_rois_fig_tsv_sbatch.py](intertask/make_rois_fig_tsv.py) 
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_rois_fig.py](intertask/make_intertask_rois_fig.py) or [make_intertask_rois_fig.sh](intertask/make_intertask_rois_fig.sh)
- [x] Make figure-specific TSV with GLM results for active vertex [make_active_vert_fig_tsv.py](intertask/make_active_vert_fig_tsv.py) 
- [x] Make figure of GLM results active vertex [make_active_vert_fig.py](intertask/make_active_vert_fig.py)
- [x] Make final statistical maps maps with pycortex [pycortex_maps_intertask.py](intertask/pycortex_maps_intertask.py) or [pycortex_maps_intertask.sh](intertask/pycortex_maps_intertask.sh)

### Eyetracking preprocessing 
- [ ] Eyetrack preprocessing [eyetrack_preproc.py](eyetracking/eyetrack_preproc.py)

### Eyetracking postprocessing
- [ ] Extract trigger timestamps [extract_triggers.py](eyetracking/extract_triggers.py)
- [ ] Extract saccades [extract_saccades.py](eyetracking/extract_saccades.py)
- [ ] Create saccade model prediction timeseries [saccade_prediction.py](eyetracking/saccade_prediction.py)
- [ ] Generate predictions timeseries and extract stats [generate_prediction.py](eyetracking/generate_prediction.py)
- [ ] Create timeseries figures [timeseries_figures.py](eyetracking/timeseries_figures.py)
- [ ] Create stats figure [stats_figures.py](eyetracking/stats_figures.py)

### Resting-state analysis
- [ ] Preprocess the resting-state data and output it in the fsLR-91k resolution [fmriprep_sbatch_ica-aroma.py] (https://github.com/mszinte/pRF_analysis/blob/main/analysis_code/preproc/functional/fmriprep_sbatch_ica-aroma.py)
- [ ] Extract motion components with ICA-AROMA using fMRIPost-AROMA [fmripost_aroma_sbatch.py](https://github.com/mszinte/pRF_analysis/blob/main/analysis_code/postproc/rest/fmripost_aroma_sbatch.py)
- [ ] Post-process and denoise the data using XCP-D [xcp-d_aroma_sbatch.py](https://github.com/mszinte/pRF_analysis/blob/main/analysis_code/postproc/rest/xcp-d_aroma_sbatch.py)
- [ ] Compute seed-based (from GLM tasks conjunction) functional connectivity on the dense timeseries with connectome workbench ()
