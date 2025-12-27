# AMBLYO7T_PRF
## About
---
Re-analysis of 7T data from strabismic patients first published here: <br/>
Schluppeck, D., Arnoldussen, D., Hussain, Z., Besle, J., Francis, S. T., & McGraw, P. V. (2025). <br/>
Strabismus and amblyopia disrupt spatial perception but not the fidelity of cortical maps in human primary visual cortex. <br/>
Vision Research, 236, 108677. https://doi.org/10.1016/j.visres.2025.108677<br/>

---
## Authors (alphabetic order): 
---
Adrien Chopin, Uriel Lascombes, Paul V McGraw, Denis Schluppek, Martin Szinte<br/>

## Links
---
- [Notes](https://docs.google.com/document/d/1ejbu9eYmVgNDWe8nTR6uHSFs7sgJo0BPWAa6hCTsnn8/edit?usp=sharing)
- [Table](https://docs.google.com/spreadsheets/d/18cAzdLURP_OE7zDU5xodlxrTOTz60x75Q7eZ-kVVhlE/edit?usp=sharing)
- [raw data](https://gin.g-node.org/schluppeck/amblyopia-data-2025/)

## Data analysis
---

### BIDS
- [x] Download raw data locally
```cd ~/temp_data/amblyopia-data-2025;gin download --content```
- [x] Upload raw data to mesocentre
```rsync -avuz --exclude='.git' --rsh='ssh -p 8822' --progress ~/temp_data/amblyopia-data-2025/ mszinte@login.mesocentre.univ-amu.fr:/scratch/mszinte/data/amblyo7T_prf/sourcedata/```
- [x] Create participants.tsv [create_participants.sh](preproc/bids/create_participant.sh)
- [x] Convert to bids [bids_conversion.py](preproc/bids/bids_conversion.py)
- [x] Create manualy event_files for concatenated runs pRFRightEye pRFLeftEye with screen settings

## MRI Data analysis
---

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] Create sagittal video [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Make cut in the brains for flattening [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh)
- [x] Flatten the cut brains [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [x] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, mask (partial recording)  [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py)
- [x] Concatenation of runs  [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### PRF Gaussian fit
- [x] Create the video of each tasks and concatenated tasks from original .mat file [create_prf_videos.ipynb](postproc/prf/fit/create_prf_videos.ipynb)
- [x] Create the visual matrix design of each tasks and concatenated tasks [create_vdm_files.ipynb](postproc/prf/fit/create_vdm_files.ipynb)
- [x] Run pRF gaussian grid fit [prf_submit_gaussfit_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_gauss_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_derivatives.py](../analysis_code/postproc/prf/postfit/compute_gauss_derivatives.py)
- [x] Make maps with pycortex [pycortex_maps_gauss.py](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.py)

##### ROIs
- [x] Draw individual ROI on fsnative data using Inkscape
- [x] Make ROIS files [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py)

##### PRF CSS fit
- [x] CSS fit [prf_submit_css_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [compute_css_stats.py](../analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](../analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [ ] Make maps with pycortex [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [ ] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [ ] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) 
- [ ] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [ ] Merge all figures [merge_fig_css.py](analysis_code/postproc/prf/postfit/merge_fig_css.py)
