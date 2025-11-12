"""
-----------------------------------------------------------------------------------------
deface_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Deface anatomical images
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: session (e.g. ses-01)
sys.argv[5]: server group (e.g. 327)
sys.argv[6]: server job or not (1 = server, 0 = terminal)
sys.argv[7]: overwrite images (0 = no, 1 = yes)
sys.argv[8]: server job or not (1 = server, 0 = terminal)
-----------------------------------------------------------------------------------------
Output(s):
Defaced images
-----------------------------------------------------------------------------------------
To run: run python commands
>> cd ~/projects/[project]/analysis_code/preproc/bids/
>> python deface_sbatch.py [main directory] [project name] [subject] [session]
                            [group] [server_project] [overwrite] [server]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analaysis/analysis_code/preproc/bids/
python deface_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 327 b327 1 1
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# imports modules
import sys
import os
import pdb
deb = pdb.set_trace

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
session_in = int(sys.argv[4])
group = sys.argv[5]
server_project = sys.argv[6]
ovewrite_in = int(sys.argv[7])
server_in = int(sys.argv[8])
hour_proc = 4
nb_procs = 8
log_dir = "{}/{}/derivatives/pp_data/logs".format(main_dir, project_dir, subject)
try: os.makedirs(log_dir)
except: pass

# define SLURM cmd
slurm_cmd = """\
#!/bin/bash
#SBATCH -p skylake
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_deface_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_deface_%N_%j_%a.out
#SBATCH -J {subject}_deface\n\n""".format(
    nb_procs=nb_procs, hour_proc=hour_proc, 
    subject=subject, log_dir=log_dir,
    server_project=server_project)

# define FSL comande 
fsl_cmd = """\
export FSLDIR='{main_dir}/{project_dir}/code/fsl'
export PATH=$PATH:$FSLDIR/bin
source $FSLDIR/etc/fslconf/fsl.sh\n""".format(main_dir=main_dir, project_dir=project_dir)

# define change mode change group comande
chmod_cmd = """chmod -Rf 771 {main_dir}/{project_dir}\n""".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = """chgrp -Rf {group} {main_dir}/{project_dir}\n""".format(main_dir=main_dir, project_dir=project_dir, group=group)


# get files
session = session_in
t1w_filename = "{}/{}/{}/{}/anat/{}_{}_T1w.nii.gz".format(main_dir, project_dir,
                                                          subject, session, subject, session)
t2w_filename = "{}/{}/{}/{}/anat/{}_{}_T2w.nii.gz".format(main_dir, project_dir, 
                                                          subject, session, subject, session)

# sh folder & file
sh_folder = "{}/{}/derivatives/pp_data/jobs".format(main_dir, project_dir, subject)
try: os.makedirs(sh_folder)
except: pass
sh_file = "{}/{}_deface.sh".format(sh_folder,subject)

of = open(sh_file, 'w')
if server_in: of.write(slurm_cmd)
if ovewrite_in == 1:
    deface_cmd_T1 = "pydeface {fn} --outfile {fn} --force --verbose\n".format(fn = t1w_filename)
    deface_cmd_T2 = "pydeface {fn} --outfile {fn} --force --verbose\n".format(fn = t2w_filename)
    
else: 
    deface_cmd_T1 = "pydeface {fn} --verbose\n".format(fn = t1w_filename)
    deface_cmd_T2 = "pydeface {fn} --verbose\n".format(fn = t2w_filename)

of.write(fsl_cmd)    
of.write(deface_cmd_T1)
of.write(deface_cmd_T2)
of.write(chmod_cmd)
of.write(chgrp_cmd)
of.close()

print(sh_file)
# Run or submit jobs
if server_in:
    os.chdir(log_dir)
    print("Submitting {} to queue".format(sh_file))
    os.system("sbatch {}".format(sh_file))
else:
    os.system("sh {}".format(sh_file))