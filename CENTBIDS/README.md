# CENTBIDS
## About
---
*We here study cortical areas involved in both vision and eye movements in 20 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** and **HCP 170k or 91k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

---
## Authors (alphabetic order): 
---
Sina KLING, Jan-Patrick STELLMANN, Martin SZINTE


## Data analysis
---
### Preprocessing
- [ ] convert MP2RANGE files to T1w

### Resting-state preprocessing
- [ ] Preprocess the resting-state data and output it in the fsLR-91k resolution [fmriprep_sbatch_ica-aroma.py](https://github.com/mszinte/pRF_analysis/blob/main/analysis_code/preproc/functional/fmriprep_sbatch_ica-aroma.py)

### Resting-state postprocessing
- [ ] Extract motion components with ICA-AROMA using fMRIPost-AROMA [fmripost_aroma_sbatch.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/rest/fmripost_aroma_sbatch.py)
- [ ] Post-process and denoise the data using XCP-D [xcp-d_aroma_sbatch.py](https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/rest/xcp-d_aroma_sbatch.py)

### Resting-state analysis
- [ ] Compute seed-based (from the tasks conjunction glm results) functional connectivity on the dense timeseries with connectome workbench ()
- [ ] Compute winner-take-all results
- [ ] Compute seed-based partial correlation
