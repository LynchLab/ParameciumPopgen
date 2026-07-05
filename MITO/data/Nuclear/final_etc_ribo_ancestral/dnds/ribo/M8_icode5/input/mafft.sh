#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -q public
#SBATCH -J mafft_aln
#SBATCH -c 8
#SBATCH -t 02-00:00
#SBATCH --mem=16G
#SBATCH -o jlogs/%x.%j.out
#SBATCH -e jlogs/%x.%j.err
#SBATCH --export=NONE

set -euo pipefail
shopt -s nullglob

module purge
module load mamba/latest
source activate pyenv

ROOT="$(pwd)"                       # run from .../dnds_analysis/input
PROT_DIR="${ROOT}/prot_seq"
OUT_DIR="${ROOT}/mafft_aligned"

mkdir -p "${OUT_DIR}" "${ROOT}/jlogs"

# Align all protein FASTAs, excluding prot_seq/single_species/
for fa in "${PROT_DIR}"/OMA*.fa; do
    base="$(basename "${fa}" .fa)"  # OMA01173
    echo "Aligning ${base} -> ${OUT_DIR}/${base}.aln.fa"
    mafft --auto --thread "${SLURM_CPUS_PER_TASK}" "${fa}" > "${OUT_DIR}/${base}.aln.fa"
done

echo "✅ All MAFFT protein alignments completed -> ${OUT_DIR}"

