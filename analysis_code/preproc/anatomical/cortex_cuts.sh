# -----------------------------------------------------------------------------------------
# cortex_cuts.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Run freewiew to make cut to the cortex (to create a cortical flatten representation)
# -----------------------------------------------------------------------------------------
# Input(s):
# $1: project directory
# $2: project name (correspond to directory)
# $3: subject name (e.g. sub-01)
# $4: mesocentre login ID
# $5: hemisphere (lh or rh)
# -----------------------------------------------------------------------------------------
# Output(s):
# 3D patch for each hemisphere
# -----------------------------------------------------------------------------------------
# To run:
# 0. TO RUN LOCALLY WITH FREEWIEW INSTALLED (not on server)
# 1. cd to function
# >> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/anatomical/
# 2. run shell command
# sh preproc/cortex_cuts.sh [main directory] [project name] [subject name] [mesocentre_ID] [hemisphere]
# -----------------------------------------------------------------------------------------
# Exemple:
# sh cortex_cuts.sh /scratch/mszinte/data MotConf sub-01 mszinte lh
# sh cortex_cuts.sh /scratch/mszinte/data MotConf sub-01 mszinte rh
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (mail@martinszinte.net)
# -----------------------------------------------------------------------------------------

# rsync to desktop (faster processing)
echo "\n>> Copying the files to the desktop"
rsync -azuv --rsh='ssh -p 8822' --progress $4@login.mesocentre.univ-amu.fr:$1/$2/derivatives/fmriprep/freesurfer/$3 ~/temp_data/

# Check + edit pial surface
echo "\n>> Proceed to the cortex cuts : "
echo "\n>> https://invibe.nohost.me/bookstack/books/preprocessing/page/cutting-inflated-brains-to-obtain-a-flattened-cortical-surface"
echo "\n>> When you are done, save the patch as '$3/surf/$5.full.patch.3d'\n"

freeview -f ~/temp_data/$3/surf/$5.inflated:annot=aparc.a2009s -layout 1 -viewport 3d

# move the file to the right place
while true; do
    read -p "Do you wish to transfer the patch to the mesocentre? (y/n) " yn
    case $yn in
        [Yy]* ) echo "\n>> Uploading the $3 patch to mesocentre";\
                rsync -avuz --rsh='ssh -p 8822' ~/temp_data/$3/surf/$5.full.patch.3d $4@login.mesocentre.univ-amu.fr:$1/$2/derivatives/fmriprep/freesurfer/$3/surf/
        break;;
        [Nn]* ) echo "\n>> No uploading of the brainmasks to mesocentre";\
                exit;;
        * ) echo "Please answer yes or no.";;
    esac
done