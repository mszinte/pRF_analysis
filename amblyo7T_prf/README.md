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
- [raw data](https://gin.g-node.org/schluppeck/amblyopia-data-2025/)

## Data analysis
---

### BIDS
- [x] Download data from [schluppeck/amblyopia-data-2025.git](gin.g-node.org/schluppeck/amblyopia-data-2025.git)
- [x] Create participants.tsv [create_participants.sh](preproc/bids/create_participant.sh)
- [x] Convert to bids sub-01 [bids_sub-01.sh](preproc/bids/bids_conversion.sh)
- [x] Create T1w from PSIR sequence [psir_mougin_adapted.ipynb](preproc/anatomical/psir_mougin_adapted.ipynb)
- [x] Deface participants t1w image [deface_sbatch.py](../analysis_code/preproc/bids/deface_sbatch.py)
- [ ] Create event_files for each task with screen settings

## MRI Data analysis
---

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Create sagittal view video before manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [x] Manual edit of brain segmentation [pial_edits.sh](../analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](../analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] Create sagittal view video before after edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [ ] Make cut in the brains for flattening [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh)
- [ ] Flatten the cut brains [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [ ] fMRIprep [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](../analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [ ] High-pass, z-score, average and leave-one-out average [preproc_end_sbatch.py](../analysis_code/preproc/functional/preproc_end_sbatch.py)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### Inter-run correlations


