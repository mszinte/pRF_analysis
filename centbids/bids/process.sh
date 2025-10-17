# copy from JP scratch
rsync -avuz /scratch/jstellmann/data/CENTBIDS/sourcedata/ /scratch/mszinte/data/centbids/sourcedata/

# copy to main folder
rsync -avuz /scratch/mszinte/data/centbids/sourcedata/ /scratch/mszinte/data/centbids/

# remove useless scans
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_FA_TB1TFL.json
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_FA_TB1TFL.nii.gz
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_TB1TFL.json
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_TB1TFL.nii.gz
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_TB1TFL_run-03.json
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_TB1TFL_run-03.nii.gz
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_TB1TFL_run-04.json
rm -f /scratch/mszinte/data/centbids/sub-2100247523/ses-01/fmap/sub-2100247523_ses-01_TB1TFL_run-04.nii.gz
rm -Rfd /scratch/mszinte/data/centbids/sub-2100247523/ses-01/log/

# add beahvioral data
rsync -avuz /Users/martin/Library/CloudStorage/Dropbox/Data/Martin/Experiments/prfexp7t/data/CENT7T/sub-2100247523/ /Users/martin/disks/meso_S/data/centbids/sub-2100247523/

# rename retino to pRF
cd /scratch/mszinte/data/centbids/
find . -name "*task-retino*" -type f | while read -r file; do
    newname="${file//task-retino/task-pRF}"
    mv "$file" "$newname"
    echo "Would rename: $file -> $newname"
done

# copy _retino_events.json
rsync -avuz /scratch/mszinte/data/RetinoMaps/task-pRF_events.json /scratch/mszinte/data/centbids/task-pRF_events.json

# copy datasetdescription.json
rsync -avuz /scratch/jstellmann/data/CENTBIDS/dataset_description.json /scratch/mszinte/data/centbids/dataset_description.json

# copy participant json and tsv
# Sina already made them earlier
rsync -avuz /scratch/jstellmann/data/CENTBIDS/participants.json /scratch/mszinte/data/centbids/participants.json
rsync -avuz /scratch/jstellmann/data/CENTBIDS/participants.tsv /scratch/mszinte/data/centbids/participants.tsv

# copy the t1w file from JP scratch
rsync -avuz /scratch/jstellmann/data/CENTBIDS/sub-2100247523/ses-01/anat/sub-2100247523_ses-01_T1w.nii.gz /scratch/mszinte/data/c
entbids/sub-2100247523/ses-01/anat/sub-2100247523_ses-01_T1w.nii.gz

# build the singularity
cd /scratch/mszinte/data/centbids/singularity/
singularity build ./fmriprep-25.1.0.simg docker://nipreps/fmriprep:25.1.0

# copy the templateflow folder
rsync -avuz /scratch/mszinte/data/RetinoMaps/code/singularity/fmriprep_tf/ /scratch/mszinte/data/centbids/code/singularity/fmriprep_tf/

# copy the output of freesurfer from JP to right place in MS scratch
# SK already clean the subject and made the flatmap
mkdir /scratch/mszinte/data/centbids/derivatives
mkdir /scratch/mszinte/data/centbids/derivatives/fmriprep
mkdir /scratch/mszinte/data/centbids/derivatives/fmriprep/freesurfer
rsync -avuz /scratch/jstellmann/data/CENTBIDS/derivatives/freesurfer/sub-2100247523/ /scratch/mszinte/data/centbids/derivatives/fmriprep/freesurfer/sub-2100247523/

# add fsaverage folders
rsync -avuz /scratch/mszinte/data/RetinoMaps/derivatives/fmriprep/freesurfer/fsaverage/ /scratch/mszinte/data/centbids/derivatives/fmriprep/freesurfer/fsaverage/
rsync -avuz /scratch/mszinte/data/RetinoMaps/derivatives/fmriprep/freesurfer/fsaverage6/ /scratch/mszinte/data/centbids/derivatives/fmriprep/freesurfer/fsaverage6/

<<<<<<< HEAD

=======
#after running pycortex import 
#copy sub-170k from retinomaps 
rsync -avuz /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/cortex/db/sub-170k /scratch/mszinte/data/centbids/derivatives/pp_data/cortex/db/sub-170k
>>>>>>> skling

# chmod/chgrp
chmod -Rf 771 /scratch/mszinte/data/centbids/
chgrp -Rf 327 /scratch/mszinte/data/centbids/