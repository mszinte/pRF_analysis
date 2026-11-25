## About
---


# Authors (alphabetic order): 
---


# Data stucture and BIDS
---

# Structural preprocessing 
---
### *Individual subject*
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format 
- [x] Download BIDSonym singularity from the GitHub Container Registry (the Docker Hub image `peerherholz/bidsonym` is no longer accessible: *requested access to the resource is denied*). This is the latest version of the `main` branch (> v0.0.6), not v0.0.4. Run it from `~` (proot `getcwd` error when launched from `/scratch`) and make sure the output folder exists:

```bash
cd ~
mkdir -p /scratch/ulascombes/data/CRASH/code/singularity
export APPTAINER_CACHEDIR=/scratch/ulascombes/.apptainer_cache
export APPTAINER_TMPDIR=/scratch/ulascombes/.apptainer_tmp
mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR
singularity build /scratch/ulascombes/data/CRASH/code/singularity/bidsonym-latest.simg docker://ghcr.io/peerherholz/bidsonym:latest
singularity run /scratch/ulascombes/data/CRASH/code/singularity/bidsonym-latest.simg --help
```

- [x] Deface participants t1w image [bidsonym_sbatch.py](../analysis_code/preproc/bids/bidsonym_sbatch.py)
  - Use `latest` as `[bidsonym_version]` argument (e.g. `python bidsonym_sbatch.py /scratch/ulascombes/data CRASH sub-01 327 b327 latest pydeface`)
  - In the latest version, `--brainextraction` only accepts `bet` (`nobrainer` was removed after v0.0.4, only used for QC): the script automatically uses `--brainextraction bet --bet_frac 0.5` for any version other than `v0.0.4`
- [ ] Download SIAM singularity (image ~10 GB, tmp ~32 GB: run from `~` with `APPTAINER_CACHEDIR`/`APPTAINER_TMPDIR` on scratch as above) (singularity build /scratch/ulascombes/data/CRASH/code/singularity/siam-latest.simg docker://romainvalabregue/siam)
- [ ] Head tissue segmentation with [SIAM](https://github.com/romainVala/SIAM) of the defaced T1w, outputs in derivatives/siam/{subject} [siam_sbatch.py](../analysis_code/preproc/anatomical/siam_sbatch.py) (settings in [siam-settings.yml](siam-settings.yml))
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

- [ ] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, anat [preproc_end.py](../analysis_code/preproc/functional/preproc_end.py) or ro run on server [preproc_end_sbatch.py](../analysis_code/preproc/functional/preproc_end_sbatch.py)
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