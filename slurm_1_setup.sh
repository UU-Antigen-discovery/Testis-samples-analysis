
#!/usr/bin/env bash
#SBATCH --job-name=hlai_setup
#SBATCH --output=logs/hlai_setup.out
#SBATCH --error=logs/hlai_setup.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

#USER CONFIGURATION
BLASTP_BIN="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis/ncbi-blast-2.17.0+/bin/blastp"
MAKEBLASTDB_BIN="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis/ncbi-blast-2.17.0+/bin/makeblastdb"
WORKDIR="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis"   # must contain *_HLAI.csv + testis_frames.fa

TARGET_FA="${WORKDIR}/testis_frames.fa"
BLASTDB_DIR="${WORKDIR}/blastdb"
FASTA_DIR="${WORKDIR}/fasta_queries"
BLAST_OUT_DIR="${WORKDIR}/blast_results_HLAI_100_perc"
HEADER_MAP="${WORKDIR}/testis_header_map.tsv"
RENAMED_FA="${BLASTDB_DIR}/testis_frames_renamed.fa"

mkdir -p logs "${BLASTDB_DIR}" "${FASTA_DIR}" "${BLAST_OUT_DIR}"

#Sanity checks
for bin in "${BLASTP_BIN}" "${MAKEBLASTDB_BIN}"; do
    [[ -x "${bin}" ]] || { echo "ERROR: not executable: ${bin}"; exit 1; }
done
[[ -f "${TARGET_FA}" ]] || { echo "ERROR: not found: ${TARGET_FA}"; exit 1; }

#Rename long headers in testis_frames.fa → short IDs + save map
echo "[setup] Renaming headers and writing header map..."

awk '
/^>/ {
    idx++
    orig = substr($0, 2)
    short = "seq_" idx
    print short "\t" orig
    print ">" short > "/dev/stderr"
    next
}
{ print > "/dev/stderr" }
' "${TARGET_FA}" > "${HEADER_MAP}" 2> "${RENAMED_FA}"

#Build BLAST protein database
echo "[setup] Running makeblastdb..."

"${MAKEBLASTDB_BIN}" \
    -in      "${RENAMED_FA}" \
    -dbtype  prot \
    -out     "${BLASTDB_DIR}/testis_frames" \
    -parse_seqids \
    -logfile "${BLASTDB_DIR}/makeblastdb.log"

echo "[setup] Database ready: ${BLASTDB_DIR}/testis_frames"

#Convert each *_HLAI.csv to a query FASTA
echo "[setup] Converting CSV files to FASTA..."

for csv in "${WORKDIR}"/*_HLAI.csv; do
    base="$(basename "${csv}" .csv)"
    out_fa="${FASTA_DIR}/${base}.fasta"

    awk -v base="${base}" '
    BEGIN { FS=";"; seq_n=0 }
    NR==1 {
        for (i=1; i<=NF; i++) {
            gsub(/^[ \t"]+|[ \t"]+$/, "", $i)
            if ($i == "Peptide") { col=i; break }
        }
        if (!col) { print "ERROR: Peptide column not found in " FILENAME > "/dev/stderr"; exit 1 }
        next
    }
    {
        val = $col
        gsub(/^[ \t"]+|[ \t"]+$/, "", val)
        if (val == "" || seen[val]++) next
        seq_n++
        print ">" base "_seq" seq_n
        print val
    }
    ' "${csv}" > "${out_fa}"

    n="$(grep -c '^>' "${out_fa}" || true)"
    echo "    ${base}.fasta  (${n} unique peptides)"
done

#Write the sample list used by the array job to index into
SAMPLE_LIST="${WORKDIR}/sample_list.txt"
ls "${FASTA_DIR}"/*_HLAI.fasta | xargs -n1 basename | sed 's/\.fasta$//' \
    > "${SAMPLE_LIST}"

echo "[setup] Sample list written: ${SAMPLE_LIST}"
echo "[setup] Done. $(wc -l < "${SAMPLE_LIST}") samples ready for array job."
