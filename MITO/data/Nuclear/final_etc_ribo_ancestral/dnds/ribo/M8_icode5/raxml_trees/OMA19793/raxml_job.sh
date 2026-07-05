#!/bin/bash
#SBATCH -q public
#SBATCH -p general
#SBATCH -J raxml_OMA19793
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA19793/raxml_OMA19793_%j.out
#SBATCH -e /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA19793/raxml_OMA19793_%j.err


module purge
module load raxml-ng-1.1.0-openmpi-np

cd "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA19793"

OMA="OMA19793"
THREADS="4"

# Resume-safe
RESUME=""
if ls OMA19793.raxml.* >/dev/null 2>&1; then
  RESUME="--redo"
fi

# -------- ONE-LINE RAxML --------
raxml-ng-mpi --all --msa "OMA19793.codon.aln.fna" --msa-format fasta --data-type DNA --model GTR+G --prefix "OMA19793.raxml" --threads 4 --tree rand{10} --bs-trees 100 

# Copy final tree for codeml
if [[ -s "OMA19793.raxml.bestTree" ]]; then
  cp -f "OMA19793.raxml.bestTree" "OMA19793.tree"
fi

echo "DONE: OMA19793 ($(date))"
