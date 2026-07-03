#!/usr/bin/env bash
set -euo pipefail

# ===============================
# INPUT / OUTPUT
# ===============================
ALN_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo/input_codon_alignment"
OUTROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo/raxml_trees"

mkdir -p "$OUTROOT"
shopt -s nullglob

# ===============================
# LOOP OVER ALL OMA
# ===============================
for ALN in "$ALN_DIR"/OMA*.codon.aln.fna; do
    OMA=$(basename "$ALN" .codon.aln.fna)
    ODIR="$OUTROOT/$OMA"
    mkdir -p "$ODIR"

    # copy alignment into its own folder
    cp -f "$ALN" "$ODIR/${OMA}.codon.aln.fna"

    JOB="$ODIR/raxml_job.sh"

    cat > "$JOB" <<EOF
#!/bin/bash
#SBATCH -q grp_mlynch11
#SBATCH --partition=general
#SBATCH -J raxml_${OMA}
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o raxml_${OMA}_%j.out
#SBATCH -e raxml_${OMA}_%j.err
#SBATCH --mail-user=\$USER@asu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module purge
module load raxml-ng-1.1.0-openmpi-np

cd "$ODIR"

# ===============================
# ONE-LINE RAxML-NG COMMAND
# ===============================
raxml-ng-mpi --all --msa "${OMA}.codon.aln.fna" --msa-format fasta --data-type DNA --model GTR+G --prefix "${OMA}" --threads 4 --tree rand{10} --bs-trees 100

echo "✅ RAxML-NG finished for ${OMA}"
EOF

    chmod +x "$JOB"
    sbatch "$JOB"
done

echo "🚀 Submitted RAxML-NG jobs for ALL OMA alignments"
echo "📁 Output directory: $OUTROOT"

