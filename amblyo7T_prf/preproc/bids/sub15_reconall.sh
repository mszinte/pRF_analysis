#!/bin/bash
#SBATCH -p skylake
#SBATCH -A b327
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH -e /scratch/mszinte/data/amblyo7T_prf/derivatives/fmriprep/log_outputs/sub-15_reconall_%N_%j.err
#SBATCH -o /scratch/mszinte/data/amblyo7T_prf/derivatives/fmriprep/log_outputs/sub-15_reconall_%N_%j.out
#SBATCH -J sub-15_reconall
#SBATCH --mail-type=NONE

export SUBJECTS_DIR=/scratch/mszinte/data/amblyo7T_prf/derivatives/fmriprep/freesurfer

# Wipe surface-derived files to ensure consistent v6.0.0 reconstruction
rm -rf $SUBJECTS_DIR/sub-15_ses-01/surf/*
rm -rf $SUBJECTS_DIR/sub-15_ses-01/label/*
rm -rf $SUBJECTS_DIR/sub-15_ses-01/stats/*

# Rerun surface reconstruction from scratch using existing mri/ volumes
recon-all -s sub-15_ses-01 -autorecon2 -autorecon3

chmod -Rf 771 /scratch/mszinte/data/amblyo7T_prf
chgrp -Rf 327 /scratch/mszinte/data/amblyo7T_prf