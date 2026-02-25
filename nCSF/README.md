## About
---
Experiment in which we use pink noise different filtered to display specific <br/>
range of spatial frequency and contrast levels to determine the neural Contrast <br/>
Sensitivity Function (nCSF) of visual cortical areas. Participants are instructed<br/>
to fixate and report an orientation of the noise pattern presented on every TR.

# Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

# Data stucture and BIDS
---
- [x] Run heudiconv on Xnat, export BIDS data and import them on meso (see : https://invibe.nohost.me/bookstack/books/7t-crmbm/page/export-data-from-xnat)
- [x] Rename subject and session [rename_subject_and_session.py](MRI/preproc/bids/rename_subject_and_session.py)
- [x] Copy events.tsv and matlab files from experients repo [copy_events_and_matfiles.py](MRI/preproc/bids/copy_events_and_matfiles.py)
- [x] Bids corections [bids_conversion.py](MRI/preproc/bids/bids_conversion.py)
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]

# Structural preprocessing 
---
### *Individual subject*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 
- [x] Download BIDSonym singularity (singularity build /scratch/mszinte/data/nCSF/code/singularity/bidsonym-v0.0.4.simg docker://peerherholz/bidsonym:v0.0.4) 
- [x] Deface participants t1w image [bidsonym_sbatch.py](../analysis_code/preproc/bids/bidsonym_sbatch.py)
- [x] Download fMRIprep singulatity (singularity build /scratch/mszinte/data/nCSF/code/singularity/fmriprep-25.2.3.simg docker://nipreps/fmriprep:25.2.3)
- [x] Download [template flow](https://github.com/templateflow/templateflow) using [datalad](https://www.datalad.org/#install) and pute it in /code/singularity/fmriprep_tf 

```bash
datalad install -r ///templateflow
cd templateflow
datalad get *
```

- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Create sagittal view video before manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Manual edit of brain segmentation [pial_edits.sh](../analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](../analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] Flattening of the cortex using [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py) (with or without manual cut [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh))

# Functional preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, anat [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py) or ro run on server [preproc_end_sbatch.py](../analysis_code/preproc/functional/preproc_end_sbatch.py)
- [x] Averaging across runs [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](preproc/functional/averaging.py)

# Inter-run correlations 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [x] Compute inter-run correlation [compute_run_corr_sbatch](../analysis_code/preproc/functional/compute_run_corr_sbatch.py)
- [x] Make maps with pycortex [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](../analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

### *Group-level analysis*
Analysis are run on the template of the HCP cifti format (**170k**) in which individual results are averaged.</br> 

- [x] Compute inter-run correlation for **template_avg** [compute_run_corr.py](../analysis_code/preproc/functional/compute_run_corr.py)
- [x] Make maps with pycortex for **hcp1.6mm** [pycortex_maps_run_corr.py](../analysis_code/preproc/functional/pycortex_maps_run_corr.py)

# pRF 
---
### *Subject-level analysis*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

#### PRF ROIs
- [x] Create 170k MMP rois masks [create_rois-mmp_npz.py](../analysis_code/atlas/create_rois-mmp_npz.py)
- [x] Project MMP atlas on fsnative surface using freesurfer [freesurfer_project_mmp_fsnative.py](../analysis_code/postproc/prf/postfit/freesurfer_project_mmp_fsnative.py)
- [x] Make fsnative roi mmp npz and overlays [fsnative_mmp_rois.py](../analysis_code/postproc/prf/postfit/fsnative_mmp_rois.py)
- [x] Make ROIS files [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

#### PRF CSS fit
- [x] Create the visual matrix design [vdm_builder_prf.py](postproc/prf/fit/vdm_builder_prf.py)
- [x] CSS fit [prf_submit_css_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [css_stats_sbatch.py](../analysis_code/postproc/prf/postfit/css_stats_sbatch.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](../analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make maps with pycortex [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [x] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) or [make_rois_fig_tsv.sh](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.sh)
- [x] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](../analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [x] Merge all figures [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)