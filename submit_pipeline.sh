#!/usr/bin/env bash
#submit_pipeline.sh  — submits all three SLURM jobs in dependency order
#
#Usage:  bash submit_pipeline.sh
#
#The sequence is as follows:
#slurm_1_setup.sh(1 job) --> slurm_2_blast_array.sh(N array tasks, one per sample)-->slurm_3_finalize.sh(1 job, runs after ALL array tasks)

set -euo pipefail

WORKDIR="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis"

#Set the directory from which the script are to be run
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${WORKDIR}/logs"

#Submit setup job
JOB1=$(sbatch \
    --parsable \
    --export=ALL \
    "${SCRIPT_DIR}/slurm_1_setup.sh")

echo "Submitted setup job: ${JOB1}"

#Wait for setup to finish to count samples and set the array range

echo "Waiting for setup job ${JOB1} to complete..."
while true; do
    state=$(squeue -j "${JOB1}" -h -o "%T" 2>/dev/null || echo "UNKNOWN")
    if [[ "${state}" == "" || "${state}" == "UNKNOWN" ]]; then
        # Job no longer in queue → finished (or failed)
        break
    fi
    sleep 30
done

#Check whether the jobs launched successfully
exit_code=$(sacct -j "${JOB1}" --format=JobID,ExitCode --noheader \
    | awk '$1 ~ /^[0-9]+$/ {gsub(/ /,"",$2); split($2,a,":"); print a[1]; exit}')
if [[ "${exit_code}" != "0" ]]; then
    echo "ERROR: Setup job ${JOB1} failed (ExitCode=${exit_code}). Aborting."
    exit 1
fi

SAMPLE_LIST="${WORKDIR}/sample_list.txt"
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "Setup complete. ${N_SAMPLES} samples found."

#Submit the BLAST array job
JOB2=$(sbatch \
    --parsable \
    --array="1-${N_SAMPLES}" \
    --export=ALL \
    "${SCRIPT_DIR}/slurm_2_blast_array.sh")

echo "Submitted array job:      ${JOB2}  (tasks 1–${N_SAMPLES})"

#Submit finalize job, depends on all array tasks completing successfully
JOB3=$(sbatch \
    --parsable \
    --dependency="afterok:${JOB2}" \
    --export=ALL \
    "${SCRIPT_DIR}/slurm_3_finalize.sh")


echo "Job chain submitted:"
echo "[${JOB1}] setup→[${JOB2}] blast array (×${N_SAMPLES})→[${JOB3}] finalize"
