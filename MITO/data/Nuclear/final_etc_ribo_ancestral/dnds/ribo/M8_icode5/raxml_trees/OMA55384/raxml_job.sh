#!/bin/bash
#SBATCH -q public
#SBATCH -p general
#SBATCH -J raxml_OMA55384
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA55384/raxml_OMA55384_%j.out
#SBATCH -e /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA55384/raxml_OMA55384_%j.err


module purge
module load raxml-ng-1.1.0-openmpi-np

cd "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees/OMA55384"

OMA="OMA55384"
THREADS="4"

# Resume-safe
RESUME=""
if ls OMA55384.raxml.* >/dev/null 2>&1; then
  RESUME="--redo"
fi

# -------- ONE-LINE RAxML --------
raxml-ng-mpi --all --msa "OMA55384.codon.aln.fna" --msa-format fasta --data-type DNA --model GTR+G --prefix "OMA55384.raxml" --threads 4 --tree rand{10} --bs-trees 100 

# Copy final tree for codeml
if [[ -s "OMA55384.raxml.bestTree" ]]; then
  cp -f "OMA55384.raxml.bestTree" "OMA55384.tree"
fi

echo "DONE: OMA55384 ($(date))"
