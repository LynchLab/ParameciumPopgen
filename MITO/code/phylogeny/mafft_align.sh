#!/bin/bash

#SBATCH -q grp_mlynch11

#SBATCH -J mafft_per_gene

#SBATCH -N 1

#SBATCH --mem=32G

#SBATCH --cpus-per-task=8

#SBATCH -t 2-00:00:00

#SBATCH -o mafft_per_gene_%j.out

#SBATCH -e mafft_per_gene_%j.err

#SBATCH --mail-user=$USER@asu.edu

#SBATCH --mail-type=END,FAIL



# === Paths ===

IN_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/phylogeny_with_outgroup/gene_fasta"

OUT_DIR="${IN_DIR}/aligned"

mkdir -p "$OUT_DIR"



# === Modules ===

module purge

source activate pyenv


# === Align each gene ===

for faa in "$IN_DIR"/*.faa; do

    gene=$(basename "$faa" .faa)

    echo "Aligning $gene ..."

    mafft --thread 8 --auto "$faa" > "$OUT_DIR/${gene}_aligned.faa"

done



echo "✅ All alignments saved in: $OUT_DIR"
