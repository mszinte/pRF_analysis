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
- [x] Make subject WEBGL with pycortex for **sub-170k** [pycortex_webgl_css.py](analysis_code/postproc/prf/webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](analysis_code/postproc/prf/webgl/publish_webgl.py)

##### GLM analysis
- [x] Run Glm for the differents tasks [glm_sbatch.py](analysis_code/postproc/glm/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](analysis_code/postproc/glm/postfit/compute_glm_stats.py)

### **Statistical analysis**
- [x] Make final statistique image [stats_final.py](analysis_code/postproc/stats/stats_final.py)
- [x] Average statisticales derivatives from all subjects in 170k template [170k_stats_averaging.py](analysis_code/postproc/stats/170k_stats_averaging.py) 
- [x] Make final satistiques maps maps with pycortex [pycortex_maps_stats_final.py](analysis_code/postproc/stats/pycortex_maps_stats_final.py)