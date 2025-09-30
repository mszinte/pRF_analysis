# centbids
## About
---
*We here study cortical areas involved in both vision and eye movements in 20 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** and **HCP 170k or 91k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

---
## Authors (alphabetic order): 
---
Sina KLING, Jan-Patrick STELLMANN, Martin SZINTE

# Main behavioral, structural MRI and functional MRI analysis
---

## BIDS
---
- [x] BIDS created by CRMBM
- [x] Deface participant by CRMBM
- [x] Copy files [process.sh](analysis_code/pRF_analysis/centbids/bids/process.sh)
- [ ] Generate experimental design matrix [create_design_matrix.py](analysis_code/preproc/bids/create_design_matrix.py)

## MRI Data analysis
---
### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] same steps in pRF_analysis

#### Functional preprocessing
- [x] same steps in pRF_analysis

## Data analysis
---

#### Functional postprocessing
