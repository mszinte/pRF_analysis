# About
---


# Authors (alphabetic order): 
---


# Data stucture and BIDS
---


# Structural preprocessing 
---
### *Subject-level analysis*


# Functional preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [ ] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [ ] High-pass, z-score, anat [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py)
- [ ] Averaging across runs [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](preproc/functional/averaging.py)

# Inter-run correlations 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [ ] Compute inter-run correlation [compute_run_corr_sbatch](../analysis_code/preproc/functional/compute_run_corr_sbatch.py)
- [ ] Make maps with pycortex [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](../analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**170k**) in which individual results are averaged.</br> 

- [ ] Compute inter-run correlation for **template_avg** [compute_run_corr.py](../analysis_code/preproc/functional/compute_run_corr.py)
- [ ] Make maps with pycortex for **hcp1.6mm** [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py)

# pRF 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

#### PRF Gaussian fit
- [ ] Create the visual matrix design [vdm_builder_prf.py](postproc/prf/fit/vdm_builder_prf.py)
- [ ] Run pRF gaussian fit [prf_submit_gaussfit_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_gauss_jobs.py)
- [ ] Compute pRF gaussian fit derivatives [compute_gauss_derivatives.py](../analysis_code/postproc/prf/postfit/compute_gauss_derivatives.py)
- [ ] Make pRF maps with pycortex [pycortex_maps_gauss.py](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.py) or [pycortex_maps_gauss.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.sh)

#### PRF ROIs
- [ ] Create 170k MMP rois masks [create_rois-mmp_npz.py](../analysis_code/atlas/create_rois-mmp_npz.py)
- [ ] Draw individual ROI on fsnative data using Inkscape and create masks rois mask [create_rois-drawn_npz.py](../analysis_code/postproc/prf/postfit/create_rois-drawn_npz.py)
- [ ] Project MMP atals on fsnative surface using freesurfer [freesurfer_project_mmp_fsnative.py](../analysis_code/postproc/prf/postfit/freesurfer_project_mmp_fsnative.py)
- [ ] Make fsnative roi mmp npz and overlays [fsnative_mmp_rois.py](../analysis_code/postproc/prf/postfit/fsnative_mmp_rois.py)
- [ ] Make ROIS files [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [ ] Create flatmaps of ROIs [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

#### PRF CSS fit
- [ ] CSS fit within the ROIs [prf_submit_css_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [ ] Compute CSS statistics [css_stats_sbatch.py](../analysis_code/postproc/prf/postfit/css_stats_sbatch.py)
- [ ] Compute CSS fit derivatives [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [ ] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](../analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [ ] Make maps with pycortex [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [ ] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [ ] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) or [make_rois_fig_tsv.sh](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.sh)
- [ ] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](../analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [ ] Merge all figures [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**170k**) in which individual results are averaged and on an ROI-based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

#### PRF Gaussian fit
- [ ] Compute pRF gaussian grid fit derivatives for **template_avg** [compute_gauss_derivatives.py](../analysis_code/postproc/prf/postfit/compute_gauss_derivatives.py)
- [ ] Make pRF maps with pycortex for **hcp1.6mm**  [pycortex_maps_gauss.py](../analysis_code/postproc/prf/postfit/pycortex_maps_gauss.py)

#### PRF ROIs
- [ ] Make ROIS files for **template_avg** [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [ ] Create flatmaps of ROIs for **hcp1.6mm** [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py)

#### PRF CSS fit
- [ ] Compute CSS statistics for **template_avg** [compute_css_stats.py](../analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [ ] Compute CSS fit derivatives for **template_avg** [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [ ] Compute CSS population cortical magnification (CM) for **template_avg** [compute_css_pcm.py](../analysis_code/postproc/prf/postfit/compute_css_pcm.py)
- [ ] Make maps with pycortex for **hcp1.6mm** [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [ ] Make general TSV with CSS pRF fit derivatives, statistics and CM for **hcp1.6mm** [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [ ] Make ROIs figure of CSS pRF fit derivatives, statistics and CM for **hcp1.6mm** and **group** [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py)
- [ ] Make ROIs figure of CSS pRF fit derivatives, statistics and CM for **hcp1.6mm** and **group** [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [ ] Merge all figures for **hcp1.6mm** and **group** [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)


# WEBGL
---
https://invibe.nohost.me/xxxxxx/

- [ ] Make subject WEBGL with pycortex [pycortex_webgl_css.py](webgl/pycortex_webgl_css.py)
- [ ] Edit [index.html](../analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](webgl/publish_webgl.py)


# pRF behaviour
---
- [ ] Make pRF behaviour figure TSV [make_prf_beh_fig_tsv.py](pRF_beh/make_prf_beh_fig_tsv.py)
- [ ] Make pRF behaviour figure [make_prf_beh_fig.py](pRF_beh/make_prf_beh_fig.py)

