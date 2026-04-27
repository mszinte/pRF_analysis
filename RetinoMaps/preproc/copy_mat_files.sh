#!/bin/bash

subjects=(
  sub-01 sub-02 sub-03 sub-04 sub-05 sub-06 sub-07 sub-08 sub-09
  sub-11 sub-12 sub-13 sub-14 sub-17 sub-20 sub-21 sub-22 sub-23 sub-24 sub-25
)

for sub in "${subjects[@]}"; do
  if [ "$sub" = "sub-01" ]; then
    ses="01"
  else
    ses="02"
  fi

  src=~/projects/locEMexp/data/${sub}/ses-${ses}/add
  dest=/scratch/mszinte/data/RetinoMaps/sourcedata/${sub}/ses-${ses}/

  if [ -d "$src" ]; then
    mkdir -p "$dest"
    rsync -av --ignore-existing "$src" "$dest"
  else
    echo "Missing: $src"
  fi
done