cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit

# ecc size pcm
sh make_ecc_size_pcm_tsv.sh /home/mszinte/projects amblyo7T_prf /scratch/mszinte/data 327
sh make_ecc_size_pcm_fig.sh /home/mszinte/projects amblyo7T_prf /scratch/mszinte/data 327

# ecc comparison
sh make_ecc_comp_tsv.sh /home/mszinte/projects amblyo7T_prf /scratch/mszinte/data 327
sh make_ecc_comp_fig.sh [code directory] [project name] [main directory] [group]

# copy all figures
sh copy_subject_figs_to_group.sh /scratch/mszinte/data amblyo7T_prf