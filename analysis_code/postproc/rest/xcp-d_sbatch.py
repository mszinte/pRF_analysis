"""
-----------------------------------------------------------------------------------------
xcp_d_sbatch.py
-----------------------------------------------------------------------------------------

Goal of the script:
Run xcp-d on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: server nb of hour to request (e.g 10)
sys.argv[5]: email account
sys.argv[6]: data group (e.g. 327)
sys.argv[7]: server_project
-----------------------------------------------------------------------------------------
Outputs: postprocessed resting-state fMRI data in a variety of formats and templates
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/analysis_code/postproc/rest
Basic command:
python xcp-d_sbatch.py /scratch/mszinte/data RetinoMaps sub-22 20 marco.bedini@univ-amu.fr 327 b327

-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr) based on the fmriprep_sbatch.py example
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
email = sys.argv[5]
group = sys.argv[6]
server_project = sys.argv[7]

# Define cluster/server specific parameters
cluster_name  = 'skylake'
singularity_img = "{main_dir}/{project_dir}/code/singularity/xcp_d-0.10.0.simg".format(
    main_dir=main_dir, project_dir=project_dir)
nb_procs = 32
memory_val = 100
log_dir = "{main_dir}/{project_dir}/derivatives/xcp-d/log_outputs".format(
    main_dir=main_dir, project_dir=project_dir)

# define SLURM cmd
slurm_cmd = """\
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -p {cluster_name}
#SBATCH --mail-user={email}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_xcp-d_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_xcp-d_%N_%j_%a.out
#SBATCH -J {subject}_xcp-d
#SBATCH --mail-type=BEGIN,END\n\n
""".format(server_project=server_project, nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val,
           log_dir=log_dir, email=email,
           cluster_name=cluster_name)

# define singularity cmd
singularity_cmd = "singularity run --cleanenv -B {main_dir}:/work_dir {simg} \
        {main_dir}/{project_dir}/derivatives/fmriprep/fmriprep_aroma {main_dir}/{project_dir}/derivatives/xcp-d/{subject} \
        participant --participant-label {sub_num} \
    	--mode none -t rest \
            --nprocs {nb_procs} --omp-nthreads {nb_procs:.0f} \
            --mem-gb {memory_val} -vvv --low-mem \
            --input-type fmriprep --smoothing 0 \
            --dummy-scans auto --despike y -p 36P \
            -w /work_dir/temp/ \
            --resource-monitor --write-graph \
            --debug all --fs-license-file {main_dir}/{project_dir}/code/freesurfer/license.txt \
            --stop-on-first-crash \
        	--abcc-qc n --combine-runs y \
        	--fd-thresh 0 --file-format cifti \
        	--lower-bpf 0.01 --upper-bpf 0.08 \
	        --linc-qc n --motion-filter-type notch \
	        --motion-filter-order 4 -r 50 \
	        --band-stop-min 12 --band-stop-max 18 \
	        --output-type interpolated \
	        --warp-surfaces-native2std y --min-coverage 0.5".format(main_dir=main_dir, 
            project_dir=project_dir, simg=singularity_img, sub_num=sub_num, subject=subject, nb_procs=nb_procs, memory_val=memory_val)

# define permission cmd
chmod_cmd = "\nchmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "\nchgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

# create sh folder and file
sh_fn = "{main_dir}/{project_dir}/derivatives/xcp-d/jobs/sub-{sub_num}_xcp_d.sh".format(
        main_dir=main_dir, sub_num=sub_num,
        project_dir=project_dir)

os.makedirs("{main_dir}/{project_dir}/derivatives/xcp-d/jobs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)
os.makedirs("{main_dir}/{project_dir}/derivatives/xcp-d/log_outputs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)

of = open(sh_fn, 'w')
of.write("{slurm_cmd}{singularity_cmd}{chmod_cmd}{chgrp_cmd}".format(
    slurm_cmd=slurm_cmd, singularity_cmd=singularity_cmd, 
    chmod_cmd=chmod_cmd, chgrp_cmd=chgrp_cmd))
of.close()

# Submit jobs
print("Submitting {sh_fn} to queue".format(sh_fn=sh_fn))
os.chdir(log_dir)
# os.system("sbatch {sh_fn}".format(sh_fn=sh_fn))
