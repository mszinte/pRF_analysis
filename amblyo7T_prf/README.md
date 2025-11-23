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

## Data analysis
---

### BIDS
- [x] Download data from [schluppeck/amblyopia-data-2025.git](gin.g-node.org/schluppeck/amblyopia-data-2025.git)
- [x] Convert to bids sub-01 [bids_sub-01.sh](preproc/bids/bids_sub-01.sh)
- [x] Deface participants t1w image [deface_sbatch.py](../analysis_code/preproc/bids/deface_sbatch.py) 
- [ ] Create event_files for each task with screen settings
- [ ] Create participants.tsv out of table1 in https://doi.org/10.1016/j.visres.2025.108677<br/>

## MRI Data analysis
---

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [ ] fMRIprep with anat-only option [fmriprep_sbatch.py](../analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] Create sagittal view video before manual edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [ ] Manual edit of brain segmentation [pial_edits.sh](../analysis_code/preproc/anatomical/pial_edits.sh)
- [ ] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](../analysis_code/preproc/anatomical/freesurfer_pial.py)
- [ ] Create sagittal view video before after edit [sagital_view.py](../analysis_code/preproc/anatomical/sagital_view.py)
- [ ] Make cut in the brains for flattening [cortex_cuts.sh](../analysis_code/preproc/anatomical/cortex_cuts.sh)
- [ ] Flatten the cut brains [flatten_sbatch.py](../analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing


#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### Inter-run correlations


