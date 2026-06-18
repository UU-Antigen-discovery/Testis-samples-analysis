
#!/usr/bin/env bash
#SBATCH --job-name=hlai_blast
#SBATCH --output=logs/hlai_blast_%a.out
#SBATCH --error=logs/hlai_blast_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-13

set -euo pipefail

#USER CONFIGURATION
BLASTP_BIN="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis/ncbi-blast-2.17.0+/bin/blastp"
WORKDIR="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis/"

BLASTDB_DIR="${WORKDIR}/blastdb"
FASTA_DIR="${WORKDIR}/fasta_queries"
BLAST_OUT_DIR="${WORKDIR}/blast_results_HLAI_100_perc"
HEADER_MAP="${WORKDIR}/testis_header_map.tsv"
SAMPLE_LIST="${WORKDIR}/sample_list.txt"

#Check for sample name
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")"

[[ -n "${SAMPLE}" ]] || { echo "ERROR: empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

echo "[array task ${SLURM_ARRAY_TASK_ID}] Processing sample: ${SAMPLE}"

FA="${FASTA_DIR}/${SAMPLE}.fasta"
RAW_OUT="${BLAST_OUT_DIR}/${SAMPLE}_blast_raw.tsv"
FINAL_OUT="${BLAST_OUT_DIR}/${SAMPLE}_blast.tsv"
LOG="${BLAST_OUT_DIR}/${SAMPLE}_blast.log"

[[ -f "${FA}" ]] || { echo "ERROR: FASTA not found: ${FA}"; exit 1; }

#Run blastp
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qcovs"

"${BLASTP_BIN}" \
    -query        "${FA}" \
    -db           "${BLASTDB_DIR}/testis_frames" \
    -out          "${RAW_OUT}" \
    -outfmt       "${OUTFMT}" \
    -task         blastp-short \
    -evalue       1e-3   \
    -word_size    3      \
    -num_threads  "${SLURM_CPUS_PER_TASK}" \
    2>"${LOG}"

#Set parameters for blastp including percentage identity and query coverage
#For the testis samples we need to have full query coverage and percentage identity

awk 'BEGIN{FS=OFS="\t"} {gsub(/ /,"",$14)} $3 >= 100 && $14 == 100' "${RAW_OUT}" > "${RAW_OUT}.filtered"

mv "${RAW_OUT}.filtered" "${RAW_OUT}"

echo "[array task ${SLURM_ARRAY_TASK_ID}] blastp finished. Restoring original headers..."

#Restore original long headers in sseqid column (col 2)
awk '
BEGIN { FS=OFS="\t" }
NR==FNR { map[$1]=$2; next }
{ $2 = ($2 in map) ? map[$2] : $2; print }
' "${HEADER_MAP}" "${RAW_OUT}" > "${FINAL_OUT}"

rm -f "${RAW_OUT}"

n_hits="$(wc -l < "${FINAL_OUT}" || true)"
echo "[array task ${SLURM_ARRAY_TASK_ID}] Done. ${n_hits} hits written to ${FINAL_OUT}"
