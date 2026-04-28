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
- [x] Create participants.tsv [create_participants.sh](preproc/bids/create_participants.sh)
- [x] Convert to bids [bids_conversion.py](preproc/bids/bids_conversion.py)
- [x] Convert to bids subject with missing raw data [bids_conversion_special.py](preproc/bids/bids_conversion_special.py)
- [x] Deface participants t1w image [bidsonym_sbatch.py](../analysis_code/preproc/bids/bidsonym_sbatch.py)
- [x] Create manualy event_files for concatenated runswith screen settings

## MRI Data analysis
---

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] Flattening of the cortex using [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py) (with or without manual cut [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh))

#### Functional preprocessing
- [x] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, mask (partial recording)  [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py)
- [x] Concatenation/Averaging of runs  [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](preproc/functional/averaging.py)

#### ROIs
- [x] Create 170k MMP rois masks [create_rois-mmp_npz.py](../analysis_code/atlas/create_rois-mmp_npz.py)
- [x] Project MMP atlas on fsnative surface using freesurfer [freesurfer_project_mmp_fsnative.py](../analysis_code/postproc/prf/postfit/freesurfer_project_mmp_fsnative.py)
- [x] Make fsnative roi mmp npz and overlays [fsnative_mmp_rois.py](../analysis_code/postproc/prf/postfit/fsnative_mmp_rois.py)
- [x] Make ROIS files [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### Visual design
- [x] Create the video of each tasks and concatenated tasks from original .mat file [create_prf_videos.ipynb](postproc/prf/fit/create_prf_videos.ipynb)
- [x] Create the visual matrix design of each tasks and concatenated tasks [create_vdm_files.ipynb](postproc/prf/fit/create_vdm_files.ipynb)
      
##### PRF Gaussian
- [x] Gaussian fit [prf_submit_gaussfit_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_gauss_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_derivatives.py](../analysis_code/postproc/prf/postfit/compute_gauss_derivatives.py)
- [x] Make maps with pycortex [pycortex_maps_gauss.py](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.py) or [pycortex_maps_gauss.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.sh)

##### PRF CSS fit
- [x] CSS fit [prf_submit_css_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [compute_css_stats.py](../analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](../analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make maps with pycortex [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) 
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [x] Merge all figures [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)

#### Main analysis
- [x] Make active verices per eye TSV [make_active_vert_tsv.py](postproc/prf/postfit/make_active_vert_tsv.py) or [make_active_vert_tsv.sh](postproc/prf/postfit/make_active_vert_tsv.sh)
- [x] Make active verices per eye TSV [make_active_vert_fig.py](postproc/prf/postfit/make_active_vert_fig.py) or [make_active_vert_fig.sh](postproc/prf/postfit/make_active_vert_fig.sh)
- [x] Make pRF parameters eye correlation TSV [make_corr_tsv.py](postproc/prf/postfit/make_corr_tsv.py) or [make_corr_tsv.sh](postproc/prf/postfit/make_corr_tsv.sh)
- [x] Make pRF parameters eye correlation figures [make_corr_fig.py](postproc/prf/postfit/make_corr_fig.py) or [make_corr_fig.sh](postproc/prf/postfit/make_corr_fig.sh)
- [x] Make pRF ecc vs size/pCM per eye TSV [make_ecc_size_pcm_tsv.py](postproc/prf/postfit/make_ecc_size_pcm_tsv.py) or [make_ecc_size_pcm_tsv.sh](postproc/prf/postfit/make_ecc_size_pcm_tsv.sh)
- [x] Make pRF ecc vs size/pCM per eye figures [make_ecc_size_pcm_fig.py](postproc/prf/postfit/make_ecc_size_pcm_fig.py) or [make_ecc_size_pcm_fig.sh](postproc/prf/postfit/make_ecc_size_pcm_fig.sh)
- [x] Copy subject figures to group folders [copy_subject_figs_to_group.py](postproc/prf/postfit/copy_subject_figs_to_group.py) or [copy_subject_figs_to_group.sh](postproc/prf/postfit/copy_subject_figs_to_group.sh)