# MOTCONF

## About
---
*General code for main analysis of MotConf project.</br>*

---
## Authors (alphabetic order): 
---
Sina KLING, Uriel LASCOMBES, Guillaume MASSON, Pascal MAMMASSIAN, Can OLUK, & Martin SZINTE

## Data analysis
---

### BIDS
- [x] Copy relevant data from PredictEye
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]
    </br>Note: for the webpage, use FireFox and wait for at least 30 min, even if nothing seems to happen.

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [ ] fMRIprep with anat-only option [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] create sagital view video before manual edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [ ] manual edit of brain segmentation [pial_edits.sh](analysis_code/preproc/anatomical/pial_edits.sh)
- [ ] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](analysis_code/preproc/anatomical/freesurfer_pial.py)
- [ ] create sagital view video before after edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [ ] make cut in the brains for flattening [cortex_cuts.sh](analysis_code/preproc/anatomical/cortex_cuts.sh)
- [ ] flatten the cut brains [flatten_sbatch.py](analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [ ] fMRIprep [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [ ] High-pass, z-score, average and leave-one-out average and correlations [preproc_end_sbatch.py](analysis_code/preproc/functional/preproc_end_sbatch.py)
- [ ] Compute vertex areas [compute_vertex_area.py](analysis_code/preproc/anatomical/compute_vertex_area.py)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### Inter-run correlations
- [ ] Compute inter-run correlation [compute_run_corr.py](analysis_code/preproc/functional/compute_run_corr.py)
- [ ] Make inter-run correlations maps with pycortex [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

##### PRF Gaussian fit
- [ ] Create the visual matrix design [vdm_builder.py](analysis_code/postproc/prf/vdm_builder.py)
- [ ] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [ ] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [ ] Make pRF maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py) or [pycortex_maps_gridfit.sh](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.sh)

##### PRF ROIs
- [ ] Copy sub-170 containing MMP rois from [RetinoMaps](https://github.com/mszinte/RetinoMaps) project [compute_gauss_gridfit_derivatives.py](https://github.com/mszinte/RetinoMaps/blob/main/analysis_code/atlas/create_170k_mmp_rois_mask.ipynb) and mask areas in the overaly that are not covered by data's field of view.
- [ ] Create 170k MMP rois masks [create_mmp_rois_atlas.py](analysis_code/atlas/create_mmp_rois_atlas.py)
- [ ] Make ROIS files [make_rois_img.py](analysis_code/postproc/prf/postfit/make_rois_img.py)
- [ ] Create flatmaps of ROIs [pycortex_maps_rois.py](analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

##### PRF CSS fit
- [ ] CSS fit within the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [ ] Compute CSS statistics [compute_css_stats.py](analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [ ] Compute CSS fit derivatives [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [ ] Compute CSS population cortical magnification [css_pcm_sbatch.py](analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [ ] Make CSS fit derivatives and pcm maps with pycortex [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [ ] Make subject WEBGL with pycortex [pycortex_webgl_css.py](analysis_code/postproc/prf/webgl/pycortex_webgl_css.py) or [pycortex_webgl_css.sh](analysis_code/postproc/prf/webgl/pycortex_webgl_css.sh)
- [ ] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](analysis_code/postproc/prf/webgl/publish_webgl.py)
- [ ] Make TSV with CSS fit derivatives, pcm and statistics [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [ ] Make pRF derivatives and pcm main figures and figure TSV [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [ ] Merge all css pycortex and pRF derivatives and pcm main figures [merge_fig_css.py](analysis_code/postproc/prf/postfit/merge_fig_css.py)
