"""
-----------------------------------------------------------------------------------------
fmripost_aroma_sbatch.py
-----------------------------------------------------------------------------------------

Goal of the script:
Run fmripost_aroma on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-01)
sys.argv[4]: server nb of hour to request (e.g 10)
sys.argv[5]: data group (e.g. 327)
sys.argv[6]: server_project
-----------------------------------------------------------------------------------------
Outputs: resting-state fMRI data denoised from motion artifacts with ICA-AROMA (aggressive, non-aggressive, or orthogonal) 
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/analysis_code/postproc/rest
Basic command:
python fmripost_aroma_sbatch.py /scratch/mszinte/data RetinoMaps sub-05 10 327 b327
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""

# Debug 
import ipdb
deb = ipdb.set_trace

# imports modules
import os
import sys
import ipdb
opj = os.path.join

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[-2:]
hour_proc = int(sys.argv[4])
group = sys.argv[5]
server_project = sys.argv[6]

# Define cluster/server specific parameters
cluster_name = 'skylake'
singularity_img = "{main_dir}/{project_dir}/code/singularity/fmripost-aroma_v0.0.8.simg".format(
    main_dir=main_dir, project_dir=project_dir)
nb_procs = 32
memory_val = 100000
log_dir = "{main_dir}/{project_dir}/derivatives/fmripost_aroma/log_outputs".format(
    main_dir=main_dir, project_dir=project_dir)

# define SLURM cmd
slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}mb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_fmripost_aroma_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_fmripost_aroma_%N_%j_%a.out
#SBATCH -J {subject}_fmripost_aroma\n\n
""".format(server_project=server_project, nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val,
           log_dir=log_dir, cluster_name=cluster_name)

# define singularity cmd
singularity_cmd = "singularity run --cleanenv -B {main_dir}:/work_dir {simg} /work_dir/{project_dir}/derivatives/fmriprep/fmriprep_aroma/ /work_dir/{project_dir}/derivatives/fmripost_aroma/ participant \
        --participant-label {sub_num} -t rest \
        --nprocs {nb_procs} --omp-nthreads {nb_procs:.0f} \
        --mem {memory_val} \
        --melodic-dimensionality -200 \
        -w /work_dir/temp/ \
        --resource-monitor --write-graph \
        --stop-on-first-crash \
	-vvv".format(main_dir=main_dir, 
            project_dir=project_dir, simg=singularity_img, sub_num=sub_num, subject=subject, nb_procs=nb_procs, memory_val=memory_val)
            
# define permission cmd
chmod_cmd = "\nchmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "\nchgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

# create sh folder and file
sh_fn = "{main_dir}/{project_dir}/derivatives/fmripost_aroma/jobs/sub-{sub_num}_fmripost_aroma.sh".format(
        main_dir=main_dir, sub_num=sub_num,
        project_dir=project_dir)

os.makedirs("{main_dir}/{project_dir}/derivatives/fmripost_aroma/jobs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)
os.makedirs("{main_dir}/{project_dir}/derivatives/fmripost_aroma/log_outputs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)

of = open(sh_fn, 'w')
of.write("{slurm_cmd}{singularity_cmd}{chmod_cmd}{chgrp_cmd}".format(
    slurm_cmd=slurm_cmd, singularity_cmd=singularity_cmd, 
    chmod_cmd=chmod_cmd, chgrp_cmd=chgrp_cmd))
of.close()

# Submit jobs
print("Submitting {sh_fn} to queue".format(sh_fn=sh_fn))
os.chdir(log_dir)
os.system("sbatch {sh_fn}".format(sh_fn=sh_fn))





