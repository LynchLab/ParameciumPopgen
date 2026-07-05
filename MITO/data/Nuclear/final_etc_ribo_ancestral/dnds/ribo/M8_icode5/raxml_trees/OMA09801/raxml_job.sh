#!/bin/bash
#SBATCH -q public
#SBATCH -p general
#SBATCH -J raxml_OMA09801
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA09801/raxml_OMA09801_%j.out
#SBATCH -e /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA09801/raxml_OMA09801_%j.err


module purge
module load raxml-ng-1.1.0-openmpi-np

cd "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA09801"

OMA="OMA09801"
THREADS="4"

# Resume-safe
RESUME=""
if ls OMA09801.raxml.* >/dev/null 2>&1; then
  RESUME="--redo"
fi

# -------- ONE-LINE RAxML --------
raxml-ng-mpi --all --msa "OMA09801.codon.aln.fna" --msa-format fasta --data-type DNA --model GTR+G --prefix "OMA09801.raxml" --threads 4 --tree rand{10} --bs-trees 100 

# Copy final tree for codeml
if [[ -s "OMA09801.raxml.bestTree" ]]; then
  cp -f "OMA09801.raxml.bestTree" "OMA09801.tree"
fi

echo "DONE: OMA09801 ($(date))"
