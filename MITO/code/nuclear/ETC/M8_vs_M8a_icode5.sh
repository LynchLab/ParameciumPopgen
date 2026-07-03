#!/bin/bash
#SBATCH -q public
#SBATCH --partition=public
#SBATCH -J M8_vs_M8a
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o M8_vs_M8a_%j.out
#SBATCH -e M8_vs_M8a_%j.err
#SBATCH --mail-user=$USER@asu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load mamba/latest
# ===============================
# LOAD CODEML (pick one)
# ===============================
source activate pyenv
# OR:
# source ~/.bashrc
# conda activate paml

# ===============================
# INPUTS
# ===============================
ALN_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo/input_codon_alignment"
TREE_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo/pruned_trees"

# ===============================
# OUTPUT (YOUR REQUEST)
# ===============================
OUTROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo/M8_vs_M8a_icode5"
mkdir -p "$OUTROOT"

RAW_TSV="$OUTROOT/m8_vs_m8a_lrt_raw.tsv"
echo -e "gene\tlnL_M8a\tlnL_M8\tstat_2deltaLL" > "$RAW_TSV"

# ===============================
# MAIN LOOP
# ===============================
shopt -s nullglob
for ALN in "$ALN_DIR"/OMA*.codon.aln.fna; do
  GENE=$(basename "$ALN" .codon.aln.fna)
  TREE="$TREE_DIR/${GENE}.tree"

  if [[ ! -f "$TREE" ]]; then
    echo "Skipping $GENE (missing tree: $TREE)"
    continue
  fi

  ALN_ABS="$(readlink -f "$ALN")"
  TREE_ABS="$(readlink -f "$TREE")"

  GOUT="$OUTROOT/$GENE"
  D_M8A="$GOUT/M8a"
  D_M8="$GOUT/M8"
  mkdir -p "$D_M8A" "$D_M8"

  echo "[$(date)] Running M8a for $GENE"

  # ----- M8a (null) -----
  cat > "$D_M8A/m8a.ctl" <<EOF
seqfile = $ALN_ABS
treefile = $TREE_ABS
outfile = ${GENE}_M8a.mlc

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
model = 0
NSsites = 8
icode = 5

fix_kappa = 0
kappa = 2

fix_omega = 1
omega = 1

cleandata = 1
ncatG = 10
Small_Diff = .5e-6
method = 0
EOF

  ( cd "$D_M8A" && codeml m8a.ctl > m8a.log 2>&1 ) || true


  echo "[$(date)] Running M8 for $GENE"

  # ----- M8 (alternative) -----
  cat > "$D_M8/m8.ctl" <<EOF
seqfile = $ALN_ABS
treefile = $TREE_ABS
outfile = ${GENE}_M8.mlc

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
model = 0
NSsites = 8
icode = 5

fix_kappa = 0
kappa = 2

fix_omega = 0
omega = 1.5

cleandata = 1
ncatG = 10
Small_Diff = .5e-6
method = 0
EOF

  ( cd "$D_M8" && codeml m8.ctl > m8.log 2>&1 ) || true


  # ===============================
  # EXTRACT lnL + COMPUTE LRT
  # ===============================
  LNL0="NA"
  LNL1="NA"

  [[ -f "$D_M8A/${GENE}_M8a.mlc" ]] && \
    LNL0=$(grep -m1 "lnL(ntime" "$D_M8A/${GENE}_M8a.mlc" | awk '{print $5}' || true)

  [[ -f "$D_M8/${GENE}_M8.mlc" ]] && \
    LNL1=$(grep -m1 "lnL(ntime" "$D_M8/${GENE}_M8.mlc" | awk '{print $5}' || true)

  [[ -z "${LNL0:-}" ]] && LNL0="NA"
  [[ -z "${LNL1:-}" ]] && LNL1="NA"

  if [[ "$LNL0" != "NA" && "$LNL1" != "NA" ]]; then
    STAT=$(python3 - <<PY
lnL0=float("$LNL0")
lnL1=float("$LNL1")
print(2.0*(lnL1-lnL0))
PY
)
  else
    STAT="NA"
  fi

  echo -e "${GENE}\t${LNL0}\t${LNL1}\t${STAT}" >> "$RAW_TSV"
done

# ===============================
# POST-PROCESS: mixture 1/2 * chi2(df=1)
# ===============================
python3 - <<PY
import math, csv
from pathlib import Path

root = Path("$OUTROOT")
tsv = root / "m8_vs_m8a_lrt_raw.tsv"

rows=[]
with tsv.open() as f:
    next(f)
    for line in f:
        g,l0,l1,stat = line.rstrip("\n").split("\t")
        try:
            x=float(stat)
            # p-value under 1/2*chi^2_1 mixture (Self & Liang; used for boundary omega=1)
            p_mix = 0.5 * math.erfc((x/2.0)**0.5)
            sig5 = "YES" if x >= 2.71 else "NO"
            sig1 = "YES" if x >= 5.41 else "NO"
            rows.append([g,l0,l1,f"{x:.6g}",f"{p_mix:.6g}",sig5,sig1])
        except:
            rows.append([g,l0,l1,"NA","NA","NA","NA"])

out = root / "m8_vs_m8a_lrt_results.csv"
with out.open("w", newline="") as f:
    w=csv.writer(f)
    w.writerow(["gene","lnL_M8a","lnL_M8","stat_2deltaLL","p_mixture_half_chi2_df1","sig_5pct","sig_1pct"])
    w.writerows(rows)

print(f"Wrote {out}")
PY
