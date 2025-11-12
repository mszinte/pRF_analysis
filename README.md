# Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Martin SZINTE


# Environment Set Up 
``` conda env create -f environment.yml ```
``` conda activate pRF_env             ```              

Then install pycortex: 

``` CFLAGS="-std=c99" pip install git+https://github.com/gallantlab/pycortex.git ```

Then install prfpy:  

``` git clone https://github.com/VU-Cog-Sci/prfpy.git ```

cd to folder 

``` python installer.py ```
everything should be fine, so to check that all version are correct: 
``` pip install -r requirements.txt ```

# To do
- [ ] improve deface_sbatch to look for all anatomical files
- [ ] see if deface really deface MP2RAGEME with noise files and others
- [ ] put the non defaced anatomy in the sourcedata folder


# Task-specific analysis
---
After pRF analysis each project has its own analysis. The project's read me can be found at : 

- RetinoMaps [README.md](RetinoMaps/README.md)
- Amblyo_prf [README.md](amblyo_prf/README.md)
- centbids [README.md](centbids/README.md)

# Main behavioral, structural MRI and functional MRI analysis
---

## BIDS
---
- [x] Deface participants t1w image [deface_sbatch.py](analysis_code/preproc/bids/deface_sbatch.py) 
    </br>Note: run script for each subject separately.
- [x] BIDS eyetracking conversion [run_eye2bids.py](analysis_code/preproc/bids/run_eye2bids.py)
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]
    </br>Note: for the webpage, use FireFox and wait for at least 30 min, even if nothing seems to happen.
- [x] Generate experimental design matrix [create_design_matrix.py](analysis_code/preproc/bids/create_design_matrix.py)

## MRI Data analysis
---

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Create sagittal view video before manual edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] Manual edit of brain segmentation [pial_edits.sh](analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] Create sagittal view video before after edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] Make cut in the brains for flattening [cortex_cuts.sh](analysis_code/preproc/anatomical/cortex_cuts.sh)
- [x] Flatten the cut brains [flatten_sbatch.py](analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [x] fMRIprep [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, average and leave-one-out average [preproc_end_sbatch.py](analysis_code/preproc/functional/preproc_end_sbatch.py)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### Inter-run correlations
- [x] Compute inter-run correlation [compute_run_corr_sbatch](analysis_code/preproc/functional/compute_run_corr_sbatch)
- [x] Make inter-run correlations maps with pycortex [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

##### PRF Gaussian fit
- [x] Create the visual matrix design [vdm_builder.py](analysis_code/postproc/prf/fit/vdm_builder.py)
- [x] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [x] Make pRF maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py) or [pycortex_maps_gridfit.sh](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.sh)

##### PRF ROIs
- [x] Create 170k MMP rois [create_hcp_rois.ipynb](analysis_code/atlas/dev/create_hcp_rois.ipynb) or copy sub-170k in pycortex folder from Retinomaps. 
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
or we ran an ROI-based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

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
