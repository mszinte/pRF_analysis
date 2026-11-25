## About
---


# Authors (alphabetic order): 
---


# Data stucture and BIDS
---

# Structural preprocessing 
---
### *Individual subject*
Analyses are run on individual participant session (**sub-0X_ses-0X**) surface (**fsnative**) or their projection on the HCP cifti format 
- [x] Download BIDSonym singularity from github, run from ~ with APPTAINER_CACHEDIR and APPTAINER_TMPDIR on scratch (singularity build /scratch/ulascombes/data/CRASH/code/singularity/bidsonym-latest.simg docker://ghcr.io/peerherholz/bidsonym:latest)
- [x] Deface participants t1w image with bidsonym latest version [bidsonym_sbatch.py](../analysis_code/preproc/bids/bidsonym_sbatch.py)
- [x] Download SIAM singularity (singularity build /scratch/ulascombes/data/CRASH/code/singularity/siam-latest.simg docker://romainvalabregue/siam)
- [ ] Head tissue segmentation of all t1w with [SIAM](https://github.com/romainVala/SIAM) on gpu [siam_sbatch.py](../analysis_code/preproc/anatomical/siam_sbatch.py)
- [ ] Create SIAM brain mask of each session (without dura mater and vessels) as fMRIprep derivative [siam_brainmask.py](../analysis_code/preproc/anatomical/siam_brainmask.py)
- [ ] Download fMRIprep singulatity (singularity build /scratch/ulascombes/data/CRASH/code/singularity/fmriprep-25.2.7.simg docker://nipreps/fmriprep:25.2.7)
- [ ] Download [template flow](https://github.com/templateflow/templateflow) using [datalad](https://www.datalad.org/#install) and pute it in /code/singularity/fmriprep_tf 

```bash
datalad install -r ///templateflow
cd templateflow
datalad get *
```

- [ ] fMRIprep with anat-only option using SIAM brain mask and sessionwise anatomy [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] Create sagittal view video of the segmentation [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [ ] Manual edit of brain segmentation only if needed [pial_edits.sh](../analysis_code/preproc/anatomical/pial_edits.sh) and [freesurfer_pial.py](../analysis_code/preproc/anatomical/freesurfer_pial.py)
- [ ] Flattening of the cortex using [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py) (with or without manual cut [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh))

# Functional preprocessing 
---
### *Subject-level analysis*
Analyses are run on individual participant session (**sub-0X_ses-0X**) surface (**fsnative**) or their projection on the HCP cifti format 

- [ ] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [ ] High-pass, z-score, anat [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py) or ro run on server [preproc_end_sbatch.py](../analysis_code/preproc/functional/preproc_end_sbatch.py)
- [ ] Averaging across runs [averaging_sbatch.py](../analysis_code/preproc/functional/averaging_sbatch.py) using [averaging.py](../analysis_code/preproc/functional/averaging.py)

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

#### PRF ROIs
- [ ] Create 170k MMP rois masks [create_rois-mmp_npz.py](../analysis_code/atlas/create_rois-mmp_npz.py)
- [ ] Project MMP atlas on fsnative surface using freesurfer [freesurfer_project_mmp_fsnative.py](../analysis_code/postproc/prf/postfit/freesurfer_project_mmp_fsnative.py)
- [ ] Make fsnative roi mmp npz and overlays [fsnative_mmp_rois.py](../analysis_code/postproc/prf/postfit/fsnative_mmp_rois.py)
- [ ] Make ROIS files [make_rois_img.py](../analysis_code/postproc/prf/postfit/make_rois_img.py)
- [ ] Create flatmaps of ROIs [pycortex_maps_rois.py](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

#### PRF CSS fit
- [ ] Create the visual matrix design [vdm_builder_prf.py](postproc/prf/fit/vdm_builder_prf.py)
- [ ] CSS fit [prf_submit_css_jobs.py](../analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [ ] Compute CSS statistics [css_stats_sbatch.py](../analysis_code/postproc/prf/postfit/css_stats_sbatch.py)
- [ ] Compute CSS fit derivatives [compute_css_derivatives.py](../analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [ ] Compute CSS population cortical magnification (CM) [css_pcm_sbatch.py](../analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [ ] Make maps with pycortex [pycortex_maps_css.py](../analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](../analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [ ] Make general TSV with CSS pRF fit derivatives, statistics and CM [make_tsv_css.py](../analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [ ] Make ROIs figure specific TSV with CSS pRF fit derivatives, statistics and CM [make_rois_fig_tsv.py](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py) or [make_rois_fig_tsv.sh](../analysis_code/postproc/prf/postfit/make_rois_fig_tsv.sh)
- [ ] Make ROIs figure of CSS pRF fit derivatives, statistics and CM [make_rois_fig.py](../analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](../analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [ ] Merge all figures [merge_fig_prf.py](../analysis_code/postproc/prf/postfit/merge_fig_prf.py)