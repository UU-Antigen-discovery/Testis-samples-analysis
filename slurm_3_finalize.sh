
#!/usr/bin/env bash
#SBATCH --job-name=hlai_finalize
#SBATCH --output=logs/hlai_finalize.out
#SBATCH --error=logs/hlai_finalize.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

#USER CONFIGURATION
WORKDIR="/hpc/local/Rocky8/uu_immunopeptidomics/testis_analysis/"

FASTA_DIR="${WORKDIR}/fasta_queries"
BLAST_OUT_DIR="${WORKDIR}/blast_results_HLAI_100_perc"
COMBINED_TSV="${WORKDIR}/combined_blast_results.tsv"
STATS_TSV="${WORKDIR}/mapping_stats.tsv"

#Combine all per-sample blast TSVs

echo "Combining per-sample results..."

{
    printf "sample\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tpeptide_sequence\tqcovs\n"

    for tsv in "${BLAST_OUT_DIR}"/*_blast.tsv; do
        base="$(basename "${tsv}" _blast.tsv)"
        awk -v sample="${base}" \
            'BEGIN{FS=OFS="\t"} NF>0 {print sample, $0}' \
            "${tsv}"
    done
} > "${COMBINED_TSV}"

total_hits="$(( $(wc -l < "${COMBINED_TSV}") - 1 ))"
echo "Combined TSV: ${COMBINED_TSV}  (${total_hits} hit rows)"

#Per-sample mapping statistics
echo "Computing mapping statistics..."

{
    printf "sample\ttotal_query_peptides\tmapped_to_single_orf\tmapped_to_multiple_isoforms_same_gene\tmapped_to_multiple_genes\tunmapped\n"

    for fa in "${FASTA_DIR}"/*_HLAI.fasta; do
        base="$(basename "${fa}" .fasta)"
        blast_tsv="${BLAST_OUT_DIR}/${base}_blast.tsv"
        total="$(grep -c '^>' "${fa}" || true)"

        if [[ ! -s "${blast_tsv}" ]]; then
            printf "%s\t%d\t0\t0\t%d\n" "${base}" "${total}" "${total}"
            continue
        fi

        awk -v total="${total}" -v sample="${base}" '
        BEGIN { FS="\t" }
        NF>0 {
            sseqid = $2
            hits[$1][sseqid] = 1
            #Extract gene ID (first ENSG token) from the ORF header
            gene = sseqid
            if (match(sseqid, /ENSG[0-9]+\.[0-9]+/))
                gene = substr(sseqid, RSTART, RLENGTH)
            genes[$1][gene] = 1
        }
        END {
            single=0; multi_same_gene=0; multi_diff_gene=0
            for (q in hits) {
                n_orfs=0;  for (s in hits[q])  n_orfs++
                n_genes=0; for (g in genes[q]) n_genes++
                if      (n_orfs  == 1) single++
                else if (n_genes == 1) multi_same_gene++
                else                   multi_diff_gene++
            }
            mapped   = single + multi_same_gene + multi_diff_gene
            unmapped = (total > mapped) ? total - mapped : 0
            printf "%s\t%d\t%d\t%d\t%d\t%d\n",
                sample, total, single, multi_same_gene, multi_diff_gene, unmapped
        }
        ' "${blast_tsv}"
    done
} > "${STATS_TSV}"

echo "Stats TSV: ${STATS_TSV}"

echo "Pipeline complete"
echo "Per-sample results : ${BLAST_OUT_DIR}/"
echo "Combined TSV: ${COMBINED_TSV}"
echo "Statistics TSV: ${STATS_TSV}"
