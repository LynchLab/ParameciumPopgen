#!/bin/bash
#SBATCH -q public
#SBATCH --partition=general
#SBATCH -J ancestral_dnds
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o ancestral_dnds_%j.out
#SBATCH -e ancestral_dnds_%j.err
#SBATCH --mail-user=$USER@asu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -u -o pipefail

BASE="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo"
PHY_DIR="$BASE/phylip_input"
TREE_DIR="$BASE/pruned_trees"
OUT_ROOT="$BASE/ancestral_dnds_icode5"

mkdir -p "$OUT_ROOT"

module purge
# module load paml/4.9j   # <-- uncomment/adjust if you use modules

if ! command -v codeml >/dev/null 2>&1; then
  echo "❌ codeml not found in PATH. Load PAML first."
  exit 1
fi

# collect all PHY files robustly
mapfile -t PHY_FILES < <(find "$PHY_DIR" -maxdepth 1 -type f -name 'OMA*.phy' | sort)

echo "Found ${#PHY_FILES[@]} PHYLIP files in: $PHY_DIR"
echo "Trees directory: $TREE_DIR"
echo "Output root:     $OUT_ROOT"
echo ""

n_ok=0
n_fail=0
n_skip=0

for phy in "${PHY_FILES[@]}"; do
  oma="$(basename "$phy" .phy)"
  tree="$TREE_DIR/${oma}.tree"

  if [[ ! -s "$phy" ]]; then
    echo "⚠️  SKIP $oma: empty phy ($phy)"
    ((n_skip++))
    continue
  fi
  if [[ ! -s "$tree" ]]; then
    echo "⚠️  SKIP $oma: missing/empty tree ($tree)"
    ((n_skip++))
    continue
  fi

  outdir="$OUT_ROOT/$oma"
  mkdir -p "$outdir"

  cp -f "$phy"  "$outdir/${oma}.phy"
  cp -f "$tree" "$outdir/${oma}.tree"

  cat > "$outdir/codeml.ctl" <<EOF
      seqfile = ${oma}.phy
     treefile = ${oma}.tree
      outfile = mlc_asr.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
   CodonFreq = 2

         model = 0
      NSsites = 0

        icode = 5

   fix_kappa = 0
       kappa = 2.0

   fix_omega = 0
       omega = 1

  cleandata = 0
 RateAncestor = 2
EOF

  echo "▶ Running codeml for $oma ..."
  # avoid hanging on "Press Enter to continue"
  if ( cd "$outdir" && printf '\n' | codeml codeml.ctl > codeml.log 2>&1 ); then
    echo "✅ $oma done"
    ((n_ok++))
  else
    echo "❌ $oma FAILED (see $outdir/codeml.log)"
    ((n_fail++))
  fi
done

echo ""
echo "================ SUMMARY ================"
echo "OK:   $n_ok"
echo "FAIL: $n_fail"
echo "SKIP: $n_skip"
echo "TOTAL: ${#PHY_FILES[@]}"
echo "========================================="

