"""
-----------------------------------------------------------------------------------------
fmriprep_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run fMRIprep on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: server nb of hour to request (e.g 10)
sys.argv[5]: anat only (1) or not (0)
sys.argv[6]: use of aroma (1) or not (0)
sys.argv[7]: use Use fieldmap-free distortion correction
sys.argv[8]: skip BIDS validation (1) or not (0)
sys.argv[9]: save cifti hcp format data with 170k vertices
sys.argv[10]: dof number (e.g. 12)
sys.argv[11]: email account
sys.argv[12]: data group (e.g. 327)
sys.argv[13]: project name (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
preprocessed files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/preproc/functional
2. run python command
python fmriprep_sbatch.py [main directory] [project name] [subject num]
                          [hour proc.] [anat_only_(y/n)] [aroma_(y/n)] [fmapfree_(y/n)] 
                          [skip_bids_val_(y/n)] [cifti_output_170k_(y/n)] [fsaverage(y/n)]
                          [dof] [email account] [group] [server_project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional
python fmriprep_sbatch.py /scratch/mszinte/data MotConf sub-01 30 anat_only_n aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 uriel.lascombes@etu.univ-amu.fr 327 b327
python fmriprep_sbatch.py /scratch/mszinte/data MotConf sub-01 30 anat_only_n aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_n 6 uriel.lascombes@etu.univ-amu.fr 327 b327
python fmriprep_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 30 anat_only_n aroma_y fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 martin.szinte@etu.univ-amu.fr 327 b327
python fmriprep_sbatch.py /scratch/mszinte/data MotConf sub-01 30 anat_only_n aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 uriel.lascombes@etu.univ-amu.fr 327 b327 fmriprep-25.1.0.simg
python fmriprep_sbatch.py /scratch/mszinte/data MotConf sub-01 30 anat_only_n aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_n 6 uriel.lascombes@etu.univ-amu.fr 327 b327 fmriprep-25.1.0.simg
python fmriprep_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 30 anat_only_n aroma_y fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 martin.szinte@etu.univ-amu.fr 327 b327 fmriprep-25.1.0.simg
python fmriprep_sbatch.py /scratch/mszinte/data centbids sub-2100247523 30 anat_only_y aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 sina.kling@univ-amu.fr 327 b327 fmriprep-25.1.0.simg
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
and Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Debug 
import ipdb
deb = ipdb.set_trace

# imports modules
import os
import sys
opj = os.path.join

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[-2:]
hour_proc = int(sys.argv[4])
anat = sys.argv[5]
aroma = sys.argv[6]
fmapfree = sys.argv[7]
skip_bids_val = sys.argv[8]
hcp_cifti_val = sys.argv[9]
fsaverage_val = sys.argv[10]
dof = int(sys.argv[11])
email = sys.argv[12]
group = sys.argv[13]
server_project = sys.argv[14]

# Define cluster/server specific parameters
cluster_name  = 'skylake'
singularity_dir = "{main_dir}/{project_dir}/code/singularity/fmriprep-22.1.1.simg".format(
    main_dir=main_dir, project_dir=project_dir)
nb_procs = 32
memory_val = 100
log_dir = "{main_dir}/{project_dir}/derivatives/fmriprep/log_outputs".format(
    main_dir=main_dir, project_dir=project_dir)

# special input
anat_only, use_aroma, use_fmapfree, anat_only_end, \
use_skip_bids_val, hcp_cifti, tf_export, tf_bind, fsaverage = '','','','','', '', '', '', ''

if anat == 'anat_only_y':
    anat_only = ' --anat-only'
    anat_only_end = '_anat'
    nb_procs = 8

if aroma == 'aroma_y':
    use_aroma = ' --use-aroma --aroma-melodic-dimensionality -200'
    aroma_end = '_aroma'

if fmapfree == 'fmapfree_y':
    use_fmapfree= ' --use-syn-sdc'

if skip_bids_val == 'skip_bids_val_y':
    use_skip_bids_val = ' --skip_bids_validation'

if hcp_cifti_val == 'cifti_output_170k_y':
    hcp_cifti = ' --cifti-output 170k'
    
if fsaverage_val == 'fsaverage_y':
    tf_export = 'export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow'
    tf_bind = "-B {main_dir}/{project_dir}/code/singularity/fmriprep_tf/:/opt/templateflow".format(
        main_dir=main_dir, project_dir=project_dir) 
    fsaverage = 'fsaverage'

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
#SBATCH -e {log_dir}/{subject}_fmriprep{anat_only_end}_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_fmriprep{anat_only_end}_%N_%j_%a.out
#SBATCH -J {subject}_fmriprep{anat_only_end}
#SBATCH --mail-type=BEGIN,END\n\n{tf_export}
""".format(server_project=server_project, nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, anat_only_end=anat_only_end, memory_val=memory_val,
           log_dir=log_dir, email=email, tf_export=tf_export,
           cluster_name=cluster_name)

#define singularity cmd
singularity_cmd = "singularity run --cleanenv {tf_bind} -B {main_dir}:/work_dir {simg} --fs-license-file /work_dir/{project_dir}/code/freesurfer/license.txt --fs-subjects-dir /work_dir/{project_dir}/derivatives/fmriprep/freesurfer/ /work_dir/{project_dir}/ /work_dir/{project_dir}/derivatives/fmriprep/fmriprep{aroma_end}/ participant --participant-label {sub_num} -w /work_dir/temp/ --bold2t1w-dof {dof} --output-spaces T1w fsnative {fsaverage} {hcp_cifti} --low-mem --mem-mb {memory_val}000 --nthreads {nb_procs:.0f} {anat_only}{use_aroma}{use_fmapfree}{use_skip_bids_val}".format(
        tf_bind=tf_bind, main_dir=main_dir, project_dir=project_dir,
        simg=singularity_dir, sub_num=sub_num, nb_procs=nb_procs,
        anat_only=anat_only, use_aroma=use_aroma, use_fmapfree=use_fmapfree,
        use_skip_bids_val=use_skip_bids_val, fsaverage = fsaverage,hcp_cifti=hcp_cifti, memory_val=memory_val,
        dof=dof, aroma_end=aroma_end)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "\nchgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

chmod_cmd_temp = "\nchmod -Rf 771 {main_dir}/temp".format(main_dir=main_dir)
chgrp_cmd_temp = "chgrp -Rf {group} {main_dir}/temp".format(main_dir=main_dir, group=group)

# create sh folder and file
sh_fn = "{main_dir}/{project_dir}/derivatives/fmriprep/jobs/sub-{sub_num}_fmriprep{anat_only_end}{aroma_end}.sh".format(
        main_dir=main_dir, sub_num=sub_num, aroma_end=aroma_end,
        project_dir=project_dir, anat_only_end=anat_only_end)

os.makedirs("{main_dir}/{project_dir}/derivatives/fmriprep/jobs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)
os.makedirs("{main_dir}/{project_dir}/derivatives/fmriprep/log_outputs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)

of = open(sh_fn, 'w')
of.write("{} \n{} \n{} \n{} \n{} \n{}".format(
    slurm_cmd, singularity_cmd, chmod_cmd, chgrp_cmd, 
    chmod_cmd_temp, chgrp_cmd_temp))
of.close()

# Submit jobs
print("Submitting {sh_fn} to queue".format(sh_fn=sh_fn))
os.chdir(log_dir)
os.system("sbatch {sh_fn}".format(sh_fn=sh_fn))