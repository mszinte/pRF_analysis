# -----------------------------------------------------------------------------------------
# pial_edits.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Run freeview to edit the segmentation brainmask localy before transfering it back to 
# the mesocentre
# -----------------------------------------------------------------------------------------
# Input(s):
# $1: project directory
# $2: project name (correspond to directory)
# $3: subject name (e.g. sub-01)
# $4: mesocentre login ID
# -----------------------------------------------------------------------------------------
# Output(s):
# edited brainmask.mgz and orignal brainmask_orog.mgz
# -----------------------------------------------------------------------------------------
# To run:
# 0. TO RUN LOCALLY WITH FREEWIEW INSTALLED (not on server)
# 1. cd to function
# >> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/anatomical/
# 2. run shell command
# sh pial_edits.sh [main directory] [project name] [subject name] [mesocentre_ID]
# -----------------------------------------------------------------------------------------
# Exemple:
# sh pial_edits.sh /scratch/mszinte/data/ RetinoMaps sub-01 mszinte
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
# -----------------------------------------------------------------------------------------

# rsync to desktop (faster processing)
echo "\n>> Copying the files to the desktop"
rsync -azuv --rsh='ssh -p 8822' --progress $4@login.mesocentre.univ-amu.fr:$1/$2/derivatives/fmriprep/freesurfer/$3 ~/temp_data/

# create a copy of the origninal brainmask
NEWFILE=~/temp_data/$3/mri/brainmask_orig.mgz
if [ -f "$NEWFILE" ]; then
    echo "\n>> A copy of original brainmask already exists: $NEWFILE"
else
    echo "\n>> Creating a copy of original brainmask: $NEWFILE"
    cp ~/temp_data/$3/mri/brainmask.mgz $NEWFILE
fi

# Check if T2.mgz exists
T2_FILE=~/temp_data/$3/mri/T2.mgz
if [ -f "$T2_FILE" ]; then
    T2_VOLUME="$T2_FILE"
    echo "\n>> T2 volume found, will be loaded in freeview"
else
    T2_VOLUME=""
    echo "\n>> No T2 volume found, skipping T2 in freeview"
fi

# Check + edit pial surface
echo "\n>> Edit the brain mask following https://invibe.nohost.me/bookstack/books/preprocessing/page/manual-edition-of-the-brain-segmentation"
echo ">> When you are done, save the brainmask and quit freeview"
freeview -v ~/temp_data/$3/mri/T1.mgz \
            ${T2_VOLUME:+$T2_VOLUME} \
            ~/temp_data/$3/mri/brainmask.mgz:opacity=0.5 \
            -f ~/temp_data/$3/surf/lh.white:edgecolor=yellow \
            ~/temp_data/$3/surf/lh.pial:edgecolor=red \
            ~/temp_data/$3/surf/rh.white:edgecolor=yellow \
            ~/temp_data/$3/surf/rh.pial:edgecolor=red

# move the file to the right place
while true; do
    read -p "Do you wish to transfer the edited brainmask to the mesocentre? (y/n) " yn
    case $yn in
        [Yy]* ) echo "\n>> Uploading of the brainmasks to mesocentre";\
                rsync -azuv --rsh='ssh -p 8822' ~/temp_data/$3/mri/brainmask.mgz $4@login.mesocentre.univ-amu.fr:$1/$2/derivatives/fmriprep/freesurfer/$3/mri/
                rsync -azuv --rsh='ssh -p 8822' ~/temp_data/$3/mri/brainmask_orig.mgz $4@login.mesocentre.univ-amu.fr:$1/$2/derivatives/fmriprep/freesurfer/$3/mri/
        break;;
        [Nn]* ) echo "\n>> No uploading of the brainmasks to mesocentre";\
                exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
