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
