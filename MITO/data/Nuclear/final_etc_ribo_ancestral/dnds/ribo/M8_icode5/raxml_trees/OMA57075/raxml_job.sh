#!/bin/bash
#SBATCH -q public
#SBATCH -p general
#SBATCH -J raxml_OMA57075
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA57075/raxml_OMA57075_%j.out
#SBATCH -e /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA57075/raxml_OMA57075_%j.err


module purge
module load raxml-ng-1.1.0-openmpi-np

cd "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA57075"

OMA="OMA57075"
THREADS="4"

# Resume-safe
RESUME=""
if ls OMA57075.raxml.* >/dev/null 2>&1; then
  RESUME="--redo"
fi

# -------- ONE-LINE RAxML --------
raxml-ng-mpi --all --msa "OMA57075.codon.aln.fna" --msa-format fasta --data-type DNA --model GTR+G --prefix "OMA57075.raxml" --threads 4 --tree rand{10} --bs-trees 100 

# Copy final tree for codeml
if [[ -s "OMA57075.raxml.bestTree" ]]; then
  cp -f "OMA57075.raxml.bestTree" "OMA57075.tree"
fi

echo "DONE: OMA57075 ($(date))"
