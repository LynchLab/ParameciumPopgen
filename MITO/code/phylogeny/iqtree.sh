#!/bin/bash

#SBATCH -q grp_mlynch11

#SBATCH -J iqtree_tetra_concat

#SBATCH -N 1

#SBATCH --mem=64G

#SBATCH --cpus-per-task=8

#SBATCH -t 2-00:00:00

#SBATCH -o iqtree/iqtree_%j.out

#SBATCH -e iqtree/iqtree_%j.err

#SBATCH --mail-user=$USER@asu.edu

#SBATCH --mail-type=END,FAIL



# === Setup ===

module purge

module load iq-tree-2.2.2.7-oy



# === Output directory ===

mkdir -p iqtree_ng



# === Run IQ-TREE in one line ===

iqtree2 -s final_concatenated_318_with_tetrahymena_cleaned_ng.faa -m LG+G --seqtype AA -B 1000 --alrt 1000 --threads 8 --prefix iqtree_ng/concat_tree



echo "✅ IQ-TREE completed. Results saved in: iqtree (prefix: concat_tree)"


