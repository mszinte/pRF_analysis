#!/bin/bash
# For every subject with a 91k_backup folder:
#   1. create 91k/rest/timeseries (nothing else, other scripts create the rest)
#   2. copy the denoised timeseries files into it
#   3. set the same group and permissions as the rest of pp_data (drwxrwx--x)
#
# Run it ON THE SERVER, not through the mounted disk.
# Usage: bash make_91k_from_backup.sh /scratch/mszinte/data/RetinoMaps/derivatives/pp_data

pp_dir=$1

for backup in $pp_dir/sub-*/91k_backup; do

    sub_dir=$(dirname $backup)
    new=$sub_dir/91k/rest/timeseries
    echo "$sub_dir"

    # 1. only the timeseries folder
    mkdir -p $new

    # 2. copy files, keeping their permissions and dates (but not the owner)
    cp --preserve=mode,timestamps $backup/rest/timeseries/* $new/

    # 3. same group as the subject folder, folders in drwxrwx--x like the others
    chgrp -R --reference=$sub_dir $sub_dir/91k
    chmod 771 $sub_dir/91k $sub_dir/91k/rest $new

done

echo "Done"