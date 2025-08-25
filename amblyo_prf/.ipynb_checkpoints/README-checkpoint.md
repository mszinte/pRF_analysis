# AMBLYO_PRF
## About
---
*We here study the organization of the cortical visual system of a population of 16 amblyopic patients and 2 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** ad **HCP 170k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

---
## Authors (alphabetic order): 
---
Adrien CHOPIN, Dennis LEVI, Uriel LASCOMBES, Jian DING, Yasha SHEYNIN, Michael SILVER, & Martin SZINTE

### Inter-group analysis
*We ran a ROI based group analysis determined individually on subject surfaces fsnative and </br>*
*by type (control (RetinoMaps) vs patient) or by amblyopia type (control vs. anisometropic/strabismic/mixed)</br>*

- [x] Make pRF derivatives and pcm main figures and figure TSV for **group-patient**, **group-patient_control**, **group-aniso**, **group-strab**, **group-mixed** [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [x] Merge all css pycortex and pRF derivatives and pcm main figures for **group-patient**, **group-patient_control**, **group-aniso**, **group-strab**, **group-mixed**  [merge_fig_css.py](analysis_code/intergroup/merge_fig_css.py)
- [x] Compute inter-group results and stats [compute_inter-group.ipynb](analysis_code/inter-group/compute_inter-group.ipynb)
- [x] Make inter-group figures [make_inter-group_fig.ipynb](analysis_code/inter-group/make_inter-group_fig.ipynb)


### Control eyetracking analysis with patients and control

- [ ] Create a visual design matrix that shift as a function of an amblyope suject eye movement dataset
- [ ] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [ ] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [ ] Make pRF maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py) or [pycortex_maps_gridfit.sh](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.sh)
- [ ] CSS fit within the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)