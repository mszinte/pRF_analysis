# Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Martin SZINTE

# Task specific analysis
---
After pRF analysis each project have is one analysis, project's read me can be find at : 

- RetinoMaps [README.md](RetinoMaps/README.md)
- Amblyo_prf [README.md](amblyo_prf/README.md)

# Main beahvioral, structural-MRI and functional-MRI analysis
---

## Behavioral analysis
---
### Individual analysis
Analyses are run on individual participant (**sub-0X**) </br>

#### Eyetracking preprocessing
- [x] BIDS conversion [run_eye2bids.py](analysis_code/preproc/bids/run_eye2bids.py)
- [ ] Extract eyetraces and Preprocessing[eyetrack_preproc.py](https://github.com/mszinte/pRF_analysis/blob/skling/RetinoMaps/eyetracking/dev/eyetrack_preproc.py)

## MRI Data analysis
---
### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] Copy relevant data from PredictEye
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/]
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] create sagital view video before manual edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] manual edit of brain segmentation [pial_edits.sh](analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] create sagital view video before after edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] make cut in the brains for flattening [cortex_cuts.sh](analysis_code/preproc/anatomical/cortex_cuts.sh)
- [x] flatten the cut brains [flatten_sbatch.py](analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [x] fMRIprep [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, average and leave-one-out average [preproc_end_sbatch.py](analysis_code/preproc/functional/preproc_end_sbatch.py)
- [x] Compute vertex areas [compute_vertex_area.py](analysis_code/preproc/anatomical/compute_vertex_area.py)

#### Eyetracking preprocessing
- [ ] BIDS conversion [run_eye2bids.py](analysis_code/preproc/bids/run_eye2bids.py)
- [ ] Generate experimental design matrix [create_design_matrix.py] (https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/create_design_matrix.ipynb)
- [ ] Extract eyetraces and Preprocessing[eyetrack_preproc.py](https://github.com/mszinte/pRF_analysis/blob/skling/RetinoMaps/eyetracking/dev/eyetrack_preproc.py)
- [ ] Extract trigger timestamps [extract_triggers.py] (https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/extract_triggers.py)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Eyetracking postprocessing
- [ ] Saccade analysis [extract_saccades.py] (https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/extract_saccades.py)
- [ ] Create experimental individual figures [generate_individual_figures.py] (https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/generate_individual_figures.py)
- [ ] Quality check [prediction_purloc.py] (https://github.com/mszinte/pRF_analysis/blob/main/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/prediction_purloc.py)

##### Inter-run correlations
- [x] Compute inter-run correlation [compute_run_corr_sbatch](analysis_code/preproc/functional/compute_run_corr_sbatch)
- [x] Make inter-run correlations maps with pycortex [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

##### PRF Gaussian fit
- [x] Create the visual matrix design [vdm_builder.py](analysis_code/postproc/prf/it/vdm_builder.py)
- [x] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [x] Make pRF maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py) or [pycortex_maps_gridfit.sh](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.sh)

##### PRF ROIs
- [x] Copy sub-170 containing MMP rois from [RetinoMaps](https://github.com/mszinte/RetinoMaps) project [compute_gauss_gridfit_derivatives.py](https://github.com/mszinte/RetinoMaps/blob/main/analysis_code/atlas/create_170k_mmp_rois_mask.ipynb) and mask areas in the overaly that are not covered by data's field of view.
- [x] Create 170k MMP rois masks [create_mmp_rois_atlas.py](analysis_code/atlas/create_mmp_rois_atlas.py)
- [x] Draw individual ROI on fsnative data using Inkscape
- [x] Make ROIS files [make_rois_img.py](analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

##### PRF CSS fit
- [x] CSS fit within the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [css_stats_sbatch.py](analysis_code/postproc/prf/postfit/css_stats_sbatch.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make CSS pRF fit derivatives and CM maps with pycortex [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [x] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) or [make_rois_fig_tsv.sh](analysis_code/postproc/prf/postfit/make_rois_fig_tsv.sh)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [x] Merge all figures [merge_fig_css.py](analysis_code/postproc/prf/postfit/merge_fig_css.py)

### Group analysis
We run either analysis on the template of the HCP cifti format (**sub-170k**) in which individual results are averaged on a template </br>
or we ran a ROI based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

#### Structural preprocessing
- [x] Compute vertex areas for **sub-170k** [compute_vertex_area.py](analysis_code/preproc/anatomical/compute_vertex_area.py)

#### Functional postprocessing

##### Inter-run correlations
- [x] Compute inter-run correlation for **sub-170k** [compute_run_corr.py](analysis_code/preproc/functional/compute_run_corr.py)
- [x] Make inter-run correlations maps with pycortex for **sub-170k** [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py)

##### PRF Gaussian fit
- [x] Compute pRF gaussian grid fit derivatives for **sub-170k** [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [x] Make pRF maps with pycortex for **sub-170k**  [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py)

##### PRF ROIs
- [x] Make ROIS files for **sub-170k** [make_rois_img.py](analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs for **sub-170k** [pycortex_maps_rois.py](analysis_code/postproc/prf/postfit/pycortex_maps_rois.py)

##### PRF CSS fit
- [x] Compute CSS statistics for **sub-170k** [compute_css_stats.py](analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives for **sub-170k** [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) for **sub-170k** [compute_css_pcm.py](analysis_code/postproc/prf/postfit/compute_css_pcm.py)
- [x] Make CSS pRF fit derivatives and CM maps with pycortex for **sub-170k** [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [x] Make general TSV with CSS pRF fit derivatives, statistics and CM for **sub-170k** [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM for **sub-170k** and **group** [make_rois_fig_tsv.py](analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM for **sub-170k** and **group** [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [x] Merge all figures for **sub-170k** and **group** [merge_fig_css.py](analysis_code/postproc/prf/postfit/merge_fig_css.py)