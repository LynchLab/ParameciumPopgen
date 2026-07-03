#!/bin/bash
#SBATCH --partition=public
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=yeast_pburs_blastp
#SBATCH --output=yeast_pburs_blastp-%J.out
#SBATCH --error=yeast_pburs_blastp-%J.err

set -euo pipefail
shopt -s nullglob

module purge
module load blast-plus-2.12.0-6z

# Run from: .../burs_blast
WD="$(pwd)"

PBURS_FASTA="${WD}/pbursaria_Dd1_mac_JV3_6_annotation_v1.0.protein.fa"
YEAST_LSU="${WD}/yeast_lsu.fasta"
YEAST_SSU="${WD}/yeast_ssu.fasta"

OUT_BLASTP="${WD}/blastp"
OUT_REV="${WD}/rev_blastp"
DBDIR="${WD}/db"

mkdir -p "${OUT_BLASTP}" "${OUT_REV}" "${DBDIR}"

THREADS="${SLURM_CPUS_PER_TASK}"
EVALUE="1e-2"

# Added qlen and slen for accurate percent calculations
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# -------------------------
# Build databases
# -------------------------

DB_PBURS="${DBDIR}/pbursaria_prot"
[[ -f "${DB_PBURS}.pin" ]] || makeblastdb -in "${PBURS_FASTA}" -dbtype prot -parse_seqids -out "${DB_PBURS}"

# Combined yeast DB (LSU + SSU)
YEAST_COMBINED="${DBDIR}/yeast_lsu_ssu_combined.fasta"
cat "${YEAST_LSU}" "${YEAST_SSU}" > "${YEAST_COMBINED}"

DB_YEAST="${DBDIR}/yeast_lsu_ssu"
[[ -f "${DB_YEAST}.pin" ]] || makeblastdb -in "${YEAST_COMBINED}" -dbtype prot -parse_seqids -out "${DB_YEAST}"

# -------------------------
# Forward BLASTP:
# Yeast (query) → pbursaria DB
# -------------------------

blastp -query "${YEAST_LSU}" -db "${DB_PBURS}" -evalue "${EVALUE}" -num_threads "${THREADS}" \
  -outfmt "${OUTFMT}" -out "${OUT_BLASTP}/yeast_lsu_vs_pbursaria.blastp.tab"

blastp -query "${YEAST_SSU}" -db "${DB_PBURS}" -evalue "${EVALUE}" -num_threads "${THREADS}" \
  -outfmt "${OUTFMT}" -out "${OUT_BLASTP}/yeast_ssu_vs_pbursaria.blastp.tab"

# -------------------------
# Reverse BLASTP:
# pbursaria (query) → yeast DB
# -------------------------

blastp -query "${PBURS_FASTA}" -db "${DB_YEAST}" -evalue "${EVALUE}" -num_threads "${THREADS}" \
  -outfmt "${OUTFMT}" -out "${OUT_REV}/pbursaria_vs_yeast_lsu_ssu.blastp.tab"

# -------------------------
# Post-process: sort by 2nd column (sseqid), add header + pct columns
# -------------------------

add_pct_and_sort() {
  local infile="$1"
  local outfile="$2"

  awk 'BEGIN {
      OFS="\t";
      print "qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","aln_over_qlen_pct","aln_over_slen_pct"
  }
  {
      qlen = $13 + 0
      slen = $14 + 0
      aln_q = (qlen > 0 ? ($4 / qlen) * 100 : 0)
      aln_s = (slen > 0 ? ($4 / slen) * 100 : 0)
      printf "%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2e\t%.1f\t%.2f\t%.2f\n",
             $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,aln_q,aln_s
  }' "${infile}" | sort -k2,2 > "${outfile}"
}

add_pct_and_sort "${OUT_REV}/pbursaria_vs_yeast_lsu_ssu.blastp.tab" \
  "${OUT_REV}/pbursaria_vs_yeast_lsu_ssu.blastp.sorted_with_pct.tsv"

echo "✅ Forward BLASTP outputs in: ${OUT_BLASTP}"
echo "✅ Reverse BLASTP outputs in: ${OUT_REV}"
echo "✅ Processed + sorted output: ${OUT_REV}/pbursaria_vs_yeast_lsu_ssu.blastp.sorted_with_pct.tsv"

