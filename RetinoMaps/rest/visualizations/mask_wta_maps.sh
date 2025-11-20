#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

ATLAS_PATH="/media/marc_be/marc_be_vault1/RetinoMaps/derivatives/func_connectivity/Retinomaps_ROIs/atlas/macro_regions"

wb_command -metric-mask "group_winner_take_all_lh.shape.gii" \
"$ATLAS_PATH/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_left_hemi.shape.gii" \
"masked/group_winner_take_all_lh_masked.shape.gii" \
    
wb_command -metric-mask "group_winner_take_all_rh.shape.gii" \
"$ATLAS_PATH/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_right_hemi.shape.gii" \
"masked/group_winner_take_all_rh_masked.shape.gii" \
