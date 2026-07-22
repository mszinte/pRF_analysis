"""
-----------------------------------------------------------------------------------------
submit_nilearn_compute_partial_corr_job.py
-----------------------------------------------------------------------------------------
Goal of the script:
Submit the nilearn partial correlation script to SLURM queue
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: group (e.g. 327)
sys.argv[4]: server project (e.g. b327)
sys.argv[5]: script name to run ("task-free" or "task-constrained")
sys.argv[6]: runtime (e.g. 1:00:00)
sys.argv[7]: covariance estimator — one of "raw", "ledoit-wolf", "graphical-lasso"
             (required for BOTH modes; task-free conditions on all 106
             bilateral MMP parcels, at least as collinear a regressor set
             as the 24 macro-regions used in task-constrained, so the same
             estimator choice applies to both)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Estimator options (both modes)
  raw             : unregularized empirical covariance (Nilearn default)
  ledoit-wolf     : Ledoit-Wolf shrinkage — analytic shrinkage intensity,
                     no cross-validation, recommended default given the
                     collinearity diagnostic run on this dataset
  graphical-lasso : GraphicalLassoCV — L1-regularized precision matrix,
                     sparsity penalty selected by cross-validation; slower,
                     produces a sparse precision matrix (a different
                     scientific claim than shrinkage — see Peterson et al.
                     2025, Imaging Neuroscience)

  Output files are tagged with the estimator name, so submitting jobs with
  different estimators never overwrites another variant's results — each
  can be run and compared independently.
-----------------------------------------------------------------------------------------
Example:
conda activate pRF_env
cd projects/pRF_analysis/RetinoMaps/rest/correlations/
python submit_nilearn_compute_partial_corr_job.py /scratch/mszinte/data RetinoMaps 327 b327 task-constrained 0:30:00 ledoit-wolf
python submit_nilearn_compute_partial_corr_job.py /scratch/mszinte/data RetinoMaps 327 b327 task-free 1:00:00 ledoit-wolf
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]
server_project = sys.argv[4]
mode = sys.argv[5]   # "task-free" or "task-constrained"
proc_time = sys.argv[6]

# Covariance estimator — now required for BOTH modes, for consistency.
# The task-free script conditions on all 106 bilateral MMP parcels (larger
# and likely at least as collinear a regressor set as the 24 macro-regions
# used in task-constrained), so the same estimator choice applies there too.
# The job script and log filenames are tagged with it so concurrent
# submissions of different estimators never collide.
VALID_ESTIMATORS = ("raw", "ledoit-wolf", "graphical-lasso")

if len(sys.argv) < 8:
    raise ValueError(
        "A covariance estimator is required as argument 7 for both modes. "
        "Accepted: {}. Example:\n"
        "  python submit_nilearn_compute_partial_corr_job.py "
        "<main_dir> <project_dir> <group> <server> <mode> "
        "<runtime> <estimator>".format(", ".join(VALID_ESTIMATORS))
    )
estimator = sys.argv[7]
if estimator not in VALID_ESTIMATORS:
    raise ValueError(
        "Invalid estimator '{}'. Accepted: {}".format(
            estimator, ", ".join(VALID_ESTIMATORS)
        )
    )

# Define which way you want to run correlations regarding the target regions
if mode == "task-free":
    script_name = "nilearn_partial_corr_seed-task_by_mmp-parcel_by_hemi.py"
elif mode == "task-constrained":
    script_name = "nilearn_partial_corr_seed-task_by_macror-task_by_hemi.py"
else:
    raise ValueError(
        "Invalid mode. Use either 'task-free' or 'task-constrained'. "
        "Got: {}".format(mode)
    )

memory_val = 50 # GB
nb_procs = 16 # number of CPUs
cluster_name = 'skylake'

# Define permission commands
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Make directories
logs_dir = '{}/{}/derivatives/pp_data/log_outputs'.format(main_dir, project_dir)
os.makedirs(logs_dir, exist_ok=True)

job_dir = '{}/{}/derivatives/pp_data/jobs'.format(main_dir, project_dir)
os.makedirs(job_dir, exist_ok=True)

# Job/log name suffix — tags the .sh filename and SLURM log filenames by
# estimator so submissions of different estimators (or task-free, which
# has no estimator) never collide on disk.
job_suffix = "_{}".format(estimator)

# SLURM header
slurm_cmd = """#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={proc_time}
#SBATCH -e {logs_dir}/job_%N_%j_%a{job_suffix}.err
#SBATCH -o {logs_dir}/job_%N_%j_%a{job_suffix}.out
""".format(
    cluster_name=cluster_name,
    server_project=server_project,
    memory_val=memory_val,
    nb_procs=nb_procs,
    proc_time=proc_time,
    logs_dir=logs_dir,
    job_suffix=job_suffix,
)

# Assuming scripts are in the same directory as this Python script
script_dir = os.path.dirname(os.path.abspath(__file__))
your_script_path = os.path.join(script_dir, script_name)

# Estimator is always passed as a positional CLI argument to the compute
# script now that both task-free and task-constrained scripts require it.
main_cmd = "python3 {} {}".format(your_script_path, estimator)

chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Create job script — tagged by estimator so concurrent submissions with
# different estimators use separate .sh files rather than overwriting
# each other before sbatch has a chance to read them.
sh_fn = "{}/all_nilearn_partial_corr{}.sh".format(job_dir, job_suffix)
with open(sh_fn, 'w') as of:
    of.write("{}\n{}\n{}\n{}\n".format(slurm_cmd, main_cmd, chmod_cmd, chgrp_cmd))

# Submit job
print("Submitting ({}) to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))