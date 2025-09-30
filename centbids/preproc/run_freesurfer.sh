#!/bin/bash

# =========================================================================
# Initial Freesurfer run using settings.json
# =========================================================================

SETTINGS_JSON="settings.json"

# Extract global variables from JSON
PROJECT_DIR=$(jq -r '.project_name' "${SETTINGS_JSON}")
FREESURFER_SINGULARITY=$(jq -r '.freesurfer_singularity' "${SETTINGS_JSON}")
SUBJECTS=$(jq -r '.subjects[]' "${SETTINGS_JSON}")
ANAT_SESSIONS=$(jq -r '.anat_session[]' "${SETTINGS_JSON}")

# Default SLURM parameters
NB_PROCS=8
HOURS=6
MEMORY="100gb"
PARTITION=$(jq -r '.cluster_name' "${SETTINGS_JSON}")
MAIL=$(jq -r '.account' "${SETTINGS_JSON}")
PROJECT_NAME=$(jq -r '.project_name' "${SETTINGS_JSON}")

# Subjects directory
SUBJECTS_DIR="${PROJECT_DIR}/derivatives/freesurfer8"
mkdir -p "${SUBJECTS_DIR}"

# Loop over subjects and sessions
for sub in $SUBJECTS; do
    sub_num=$(echo $sub | sed 's/sub-//')
    for ses in $ANAT_SESSIONS; do
        echo "Processing ${sub} ${ses} ..."

        # Freesurfer input T1
        T1_SRC=$(ls "${PROJECT_DIR}/derivatives/presurfer/${sub}/${ses}/presurf_MPRAGEise/${sub}_${ses}"*".nii" 2>/dev/null | head -n 1)
        if [ -z "$T1_SRC" ]; then
            echo "Warning: No MPRAGE file found for ${sub} ${ses} - skipping"
            continue
        fi

        T1_BASE=$(basename "$T1_SRC")
        FS_ID="${sub}_${ses}"
        FS_SUBJ_DIR="${SUBJECTS_DIR}/${FS_ID}"

        # Skip if Freesurfer output already exists
        if [ -d "$FS_SUBJ_DIR" ]; then
            echo "Freesurfer already run for ${FS_ID} - skipping"
            continue
        fi

        LOG_DIR="${PROJECT_DIR}/${sub}/${ses}/log"
        mkdir -p "$LOG_DIR"
        LOGFILE="${LOG_DIR}/freesurfer_${FS_ID}.log"
        echo "Logfile: $LOGFILE"

        # Freesurfer singularity command
        SINGULARITY_CMD="apptainer exec \
--bind ${PROJECT_DIR}/derivatives/presurfer/${sub}/${ses}/presurf_MPRAGEise/:/ses/,${SUBJECTS_DIR}:\$HOME/freesurfer-subjects-dir \
--env FS_ALLOW_DEEP=1 ${FREESURFER_SINGULARITY} \
recon-all -i /ses/${T1_BASE} -s ${FS_ID} -all"

        # Prepare SLURM job file
        JOBFILE="${LOG_DIR}/freesurfer_${FS_ID}_slurm.sh"
        {
            echo "#!/bin/bash"
            echo "#SBATCH -J freesurfer_${FS_ID}"
            echo "#SBATCH --mail-type=ALL"
            echo "#SBATCH --mail-user=${MAIL}"
            echo "#SBATCH -p ${PARTITION}"
            echo "#SBATCH -A ${PROJECT_NAME}"
            echo "#SBATCH --nodes=1"
            echo "#SBATCH --cpus-per-task=${NB_PROCS}"
            echo "#SBATCH --mem=${MEMORY}"
            echo "#SBATCH --time=${HOURS}:00:00"
            echo "#SBATCH -e ${LOGFILE%.log}.err"
            echo "#SBATCH -o ${LOGFILE%.log}.out"
            echo "#SBATCH --mail-type=BEGIN,END"
            echo "$SINGULARITY_CMD"
        } > "$JOBFILE"

        # Submit the job
        echo "Submitting Freesurfer job for ${FS_ID}"
        sbatch "$JOBFILE"
    done
done
