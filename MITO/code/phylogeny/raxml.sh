#!/bin/bash

#SBATCH -q grp_mlynch11

#SBATCH -J raxml_all_species

#SBATCH -N 1

#SBATCH --mem=128G

#SBATCH --cpus-per-task=8

#SBATCH -t 7-00:00:00

#SBATCH -o raxml_all_species_%j.out

#SBATCH -e raxml_all_species_%j.err

#SBATCH --mail-user=$USER@asu.edu

#SBATCH --mail-type=BEGIN,END,FAIL



# === Load RAxML-NG ===

module purge

module load raxml-ng-1.1.0-openmpi-np



# === Paths ===

INFA="final_concatenated_318_with_tetrahymena_cleaned_ng.faa"

OUTDIR="raxml_tree_ng"

PREFIX="${OUTDIR}/bmge_raxml"



mkdir -p "$OUTDIR"



# === Run RAxML-NG in one line ===

raxml-ng-mpi --all --msa "$INFA" --msa-format fasta --model LG+G --data-type AA --prefix "$PREFIX" --threads 8 --tree rand{10} --bs-trees 100 --redo



echo "✅ RAxML-NG completed. Results saved in: $OUTDIR"


