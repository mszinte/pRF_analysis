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
##### Make project WEBGL
- [x] Make subject WEBGL with pycortex for **sub-170k** [pycortex_webgl_css.py](webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](webgl/publish_webgl.py)

##### GLM analysis
- [x] Run Glm for the differents tasks [glm_sbatch.py](glm/fit/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](glm/postfit/compute_glm_stats.py)

### Inter task analysis
- [x] Make intertasks image [make_intertask_img.py](intertask/make_intertask_img.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics, CM and GLM results [make_intertask_rois_fig_tsv_sbatch.py](intertask/make_intertask_rois_fig_tsv_sbatch.py)
- [] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [] Make final satistiques maps maps with pycortex [pycortex_maps_stats_final.py](analysis_code/postproc/stats/pycortex_maps_stats_final.py)