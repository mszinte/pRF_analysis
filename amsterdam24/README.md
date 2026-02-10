# About
---
*We here study population receptive fields for 3 healthy subjects during normal vision, and 3 altered vision conditions (pRFAniso: glasses to have vision less clear. pRFStrab: prism glasses to imitate strabism. pRFOccl: patch on one eye 3 hours advance, taken off just before scanning) </br>*


# Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

preprocesssing overview see: [https://docs.google.com/spreadsheets/d/1yQR9NYtAlSiqjRlXezVrVn0zPcKVd6oCFOHrvTpO4X8/edit?gid=0#gid=0]
experiment overview see: [https://docs.google.com/spreadsheets/d/1Y0ZTe_oxl2hHWhh2eh3WTmuDQO_HNvoqlqjh2_CPUmM/edit?gid=1491380407#gid=1491380407]

# Data stucture and BIDS
---
- [x] Convert to nifti with dcm2niix, copy relevant data and rename according to BIDS standards [process.sh](process.sh)
- [x] Copy events.tsv and matlab files from experients repo [copy_events_and_matfiles.py](preproc/copy_events_and_matfiles.py)
- [x] BIDS eyetracking conversion [run_eye2bids.py](../analysis_code/preproc/bids/run_eye2bids.py)
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]

# Structural preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Create sagittal view video before manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Manual edit of brain segmentation [pial_edits.sh](../analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](../analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] Create sagittal view video after manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Flattening of the cortex using [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py) (with or without manual cut [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh))

# Functional preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, anat [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py)
- [ ] Averaging across runs [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](preproc/functional/averaging.py)

# Inter-run correlations 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [ ] Compute inter-run correlation [compute_run_corr_sbatch](../analysis_code/preproc/functional/compute_run_corr_sbatch.py)
- [ ] Make maps with pycortex [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](../analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

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

