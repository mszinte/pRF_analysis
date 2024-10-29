# Copy main bids folder from PredictEye project
rsync -avuz --progress /scratch/mszinte/data/PredictEye/bids_data/ /scratch/mszinte/data/RetinoMaps

# Create additional BIDS folders
mkdir /scratch/mszinte/data/RetinoMaps/code
mkdir /scratch/mszinte/data/RetinoMaps/derivatives
mkdir /scratch/mszinte/data/RetinoMaps/sourcedata

# Copy code folder from amblyo_prf
rsync -avuz --progress /scratch/mszinte/data/amblyo_prf/code/ /scratch/mszinte/data/RetinoMaps/code

# Change permissions to all data
chmod -Rf 771 /scratch/mszinte/data/RetinoMaps
chgrp -Rf 327 /scratch/mszinte/data/RetinoMaps