#!/usr/bin/env bash

# ===============================
# PATHS
# ===============================
ALN_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/input/4_codon_alignment"
OUTROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/raxml_trees"

mkdir -p "$OUTROOT"

# ===============================
# SETTINGS
# ===============================
THREADS=4
MEM="16G"
TIME="1-00:00:00"
N_START_TREES=10
N_BOOTSTRAP=100

# ===============================
# LOOP OVER ALIGNMENTS
# ===============================
for ALN in "$ALN_DIR"/OMA*.codon.aln.fna; do
  OMA=$(basename "$ALN" .codon.aln.fna)

  # Count species
  NSPECIES=$(grep -c '^>' "$ALN")

  if [[ "$NSPECIES" -le 2 ]]; then
    echo "✗ Skip $OMA (only $NSPECIES species)"
    continue
  fi

  ODIR="$OUTROOT/$OMA"
  mkdir -p "$ODIR"

  # Skip if tree already exists
  if [[ -s "$ODIR/${OMA}.tree" || -s "$ODIR/${OMA}.raxml.bestTree" ]]; then
    echo "✓ Skip $OMA (tree exists)"
    continue
  fi

  cp -f "$ALN" "$ODIR/${OMA}.codon.aln.fna"

  JOB="$ODIR/raxml_job.sh"
  cat > "$JOB" <<EOF
#!/bin/bash
#SBATCH -q public
#SBATCH -p general
#SBATCH -J raxml_${OMA}
#SBATCH -N 1
#SBATCH --mem=${MEM}
#SBATCH --cpus-per-task=${THREADS}
#SBATCH -t ${TIME}
#SBATCH -o ${ODIR}/raxml_${OMA}_%j.out
#SBATCH -e ${ODIR}/raxml_${OMA}_%j.err


module purge
module load raxml-ng-1.1.0-openmpi-np

cd "$ODIR"

OMA="${OMA}"
THREADS="${THREADS}"

# Resume-safe
RESUME=""
if ls ${OMA}.raxml.* >/dev/null 2>&1; then
  RESUME="--redo"
fi

# -------- ONE-LINE RAxML --------
raxml-ng-mpi --all --msa "${OMA}.codon.aln.fna" --msa-format fasta --data-type DNA --model GTR+G --prefix "${OMA}.raxml" --threads ${THREADS} --tree rand{${N_START_TREES}} --bs-trees ${N_BOOTSTRAP} ${RESUME}

# Copy final tree for codeml
if [[ -s "${OMA}.raxml.bestTree" ]]; then
  cp -f "${OMA}.raxml.bestTree" "${OMA}.tree"
fi

echo "DONE: ${OMA} (\$(date))"
EOF

  chmod +x "$JOB"
  sbatch "$JOB"
  echo "→ Submitted $OMA ($NSPECIES species)"
done

echo "🚀 All eligible RAxML jobs submitted"

