# RetinoMaps
## About
---
*We here study cortical areas involved in both vision and eye movements in 20 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** ad **HCP 170k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

---
## Authors (alphabetic order): 
---
Sina KLING, Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

## To Do 
---
- [ ] take care of glm analysis

## Data analysis
---
### Prepocessing
- [x] Copy relevant data from PredictEye [copy_data.py](preproc/bids_copy_data.sh) 
- [x] Change the 'task' to 'task_condition' coulumn name in event.tsv files to avoid BIDS problems [correct_events_files.ipynb](preproc/correct_events_files.ipynb)

### pRF behaviour analysis
- [x] Make pRF behaviour figure TSV [make_prf_beh_fig_tsv.py](pRF_beh/make_prf_beh_fig_tsv.py)
- [x] Make pRF behaviour figure [make_prf_beh_fig.py](pRF_beh/make_prf_beh_fig.py)

### Make project WEBGL
- [x] Make subject WEBGL with pycortex for **sub-170k** [pycortex_webgl_css.py](webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](webgl/publish_webgl.py)

### GLM analysis
- [x] Run Glm for the differents tasks [glm_sbatch.py](glm/fit/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](glm/postfit/compute_glm_stats.py)

### Inter task analysis
- [x] Make intertasks image [make_intertask_img.py](intertask/make_intertask_img.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_tsv.py](intertask/make_intertask_tsv.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics, CM  and GLM results [make_intertask_rois_fig_tsv_sbatch.py](intertask/make_rois_fig_tsv.py) 
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_rois_fig.py](intertask/make_intertask_rois_fig.py) or [make_intertask_rois_fig.sh](intertask/make_intertask_rois_fig.sh)
- [x] Make final satistiques maps maps with pycortex [pycortex_maps_intertask.py](intertask/pycortex_maps_intertask.py) or [pycortex_maps_intertask.sh](intertask/pycortex_maps_intertask.sh)

#### Eyetracking preprocessing 
- [ ] Generate experimental design matrix [eyetracking/dev/create_design_matrix.ipynb](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/create_design_matrix.ipynb)
- [ ] Extract eyetraces and Preprocessing[eyetrack_preproc.py](eyetracking/dev/eyetrack_preproc.py)


#### Eyetracking postprocessing
- [ ] Extract trigger timestamps [eyetracking/dev/PurLoc_SacLoc/extract_triggers.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/extract_triggers.py)
- [ ] Extract saccades [eyetracking/dev/PurLoc_SacLoc/extract_saccades.py] (https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/extract_saccades.py)
- [ ] Create experimental individual figures [eyetracking/dev/PurLoc_SacLoc/generate_individual_figures.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/generate_individual_figures.py)
- [ ] Quality check PurLoc [eyetracking/dev/PurLoc_SacLoc/prediction_purloc.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/prediction_purloc.py)
- [ ] Quality check pRF [eyetracking/dev/PurLoc_SacLoc/prediction_purloc.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/prediction_pRF.py)
- [ ] Generate saccade model for Quality check SacLoc [eyetracking/dev/PurLoc_SacLoc/prediction_sacloc.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/prediction_sacloc.py)
- [ ] Quality check pRF [eyetracking/dev/PurLoc_SacLoc/prediction_sacloc.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/prediction_sacloc.py)