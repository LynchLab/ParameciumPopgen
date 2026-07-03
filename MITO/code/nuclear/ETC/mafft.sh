#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -q public
#SBATCH -J mafft_prot
#SBATCH -c 8
#SBATCH -t 02-00:00
#SBATCH --mem=16G
#SBATCH -o jlogs/%x.%j.out
#SBATCH -e jlogs/%x.%j.err
#SBATCH --export=NONE


# -------------------------
# environment
# -------------------------
module purge
module load mamba/latest
source activate pyenv

# -------------------------
# directories
# -------------------------
INPUT_DIR="protein_seq"
OUT_DIR="mafft_prot"

mkdir -p "${OUT_DIR}" jlogs

# -------------------------
# MAFFT alignment
# -------------------------
for fa in "${INPUT_DIR}"/OMA*.fa; do
    base=$(basename "${fa}" .fa)
    echo "Aligning ${base}"
    mafft --auto --thread "${SLURM_CPUS_PER_TASK}" "${fa}" > "${OUT_DIR}/${base}.aln.fa"
done

echo "✅ All MAFFT protein alignments completed"

