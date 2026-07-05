#!/bin/bash
#SBATCH -q public
#SBATCH --partition=general
#SBATCH -J burs_asr
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH -o burs_asr_%j.out
#SBATCH -e burs_asr_%j.err
#SBATCH --mail-user=$USER@asu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ===============================
# PATHS (YOUR DATA)
# ===============================
BASE="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis"

ALN_DIR="$BASE/input/4_codon_alignment"
TREE_DIR="$BASE/input/pruned_trees"

OUT_ROOT="$BASE/ancestral_dnds"
mkdir -p "$OUT_ROOT"

# Minimum taxa required to run codeml
MIN_TAXA=3

# ===============================
# ENV / CODEML
# ===============================
# NOTE: keep this if your cluster uses conda activate; otherwise replace with:
# source /path/to/miniconda3/etc/profile.d/conda.sh && conda activate pyenv
source activate pyenv 2>/dev/null || true
module purge

if ! command -v codeml >/dev/null 2>&1; then
  echo "❌ codeml not found in PATH. Load PAML first (conda/module)."
  exit 1
fi

# ===============================
# COLLECT INPUTS
# ===============================
shopt -s nullglob
ALNS=( "$ALN_DIR"/OMA*.codon.aln.fna )

echo "Alignments: $ALN_DIR  (n=${#ALNS[@]})"
echo "Trees:      $TREE_DIR"
echo "Output:     $OUT_ROOT"
echo "MIN_TAXA:   $MIN_TAXA"
echo ""

n_ok=0
n_fail=0
n_skip=0

# ===============================
# MAIN LOOP
# ===============================
for aln in "${ALNS[@]}"; do
  oma="$(basename "$aln" .codon.aln.fna)"
  tree_in="$TREE_DIR/${oma}.tree"

  if [[ ! -s "$aln" ]]; then
    echo "⚠️  SKIP $oma: empty alignment ($aln)"
    ((n_skip+=1))
    continue
  fi
  if [[ ! -s "$tree_in" ]]; then
    echo "⚠️  SKIP $oma: missing/empty tree ($tree_in)"
    ((n_skip+=1))
    continue
  fi

  outdir="$OUT_ROOT/$oma"
  mkdir -p "$outdir"

  # Build PHYLIP + clean tree that matches PHYLIP taxa
  if ! python3 - <<PY
import re, sys
from pathlib import Path

oma = "$oma"
aln = Path("$aln")
tree_in = Path("$tree_in")
outdir = Path("$outdir")

MIN_TAXA = int("$MIN_TAXA")

# ---------- Read FASTA ----------
def read_fasta(path: Path):
  records=[]
  name=None
  seq=[]
  for line in path.read_text().splitlines():
    line=line.strip()
    if not line:
      continue
    if line.startswith(">"):
      if name is not None:
        records.append((name,"".join(seq)))
      name=line[1:].split()[0]
      seq=[]
    else:
      seq.append(line)
  if name is not None:
    records.append((name,"".join(seq)))
  return records

recs = read_fasta(aln)
if not recs:
  print(f"[SKIP] {oma}: no FASTA records", file=sys.stderr)
  sys.exit(20)

# ---------- Keep only taxa with non-empty sequences ----------
recs2=[]
for n,s in recs:
  s=s.strip()
  if s and set(s) != {"-"}:
    recs2.append((n,s))
recs = recs2

if len(recs) < MIN_TAXA:
  print(f"[SKIP] {oma}: taxa={len(recs)} < MIN_TAXA={MIN_TAXA}", file=sys.stderr)
  sys.exit(21)

# ---------- Convert names to PHYLIP-friendly (<=10 chars) ----------
taxa=[]
seqs=[]
for n,s in recs:
  t = n.split("|")[0]     # PBUR|... -> PBUR  (or PBUR already)
  t = re.sub(r"[^A-Za-z0-9_.-]", "", t)
  t = t[:10]
  if not t:
    print(f"[SKIP] {oma}: bad taxon name from {n}", file=sys.stderr)
    sys.exit(22)
  taxa.append(t)
  seqs.append(s.upper())

# Ensure all same length
L = len(seqs[0])
if any(len(x)!=L for x in seqs):
  print(f"[SKIP] {oma}: sequences not equal length", file=sys.stderr)
  sys.exit(23)

# Write sequential PHYLIP
phy = outdir / f"{oma}.phy"
lines=[f"{len(taxa)} {L}"]
for t,s in zip(taxa, seqs):
  lines.append(f"{t.ljust(10)} {s}")
phy.write_text("\n".join(lines) + "\n")

# ---------- Read + sanitize tree ----------
tree = tree_in.read_text().strip()
tree = re.sub(r"\s+", "", tree)

if tree == ";" or tree == "":
  print(f"[SKIP] {oma}: tree is empty (';')", file=sys.stderr)
  sys.exit(24)

# Fix unary nodes iteratively: ((X):len) -> (X)
changed=True
while changed:
  new = re.sub(r"\(\(([^()]+)\):[0-9.eE+-]+\)", r"(\1)", tree)
  changed = (new != tree)
  tree = new

# Clean artifacts
tree = re.sub(r",\)", ")", tree)
tree = re.sub(r"\(,", "(", tree)
tree = re.sub(r",,+", ",", tree)
tree = tree.rstrip(";") + ";"

# Extract tip labels
tips = set(re.findall(r"([A-Za-z0-9_.-]+):[0-9.eE+-]+", tree))
taxa_set = set(taxa)

if tips != taxa_set:
  missing = sorted(taxa_set - tips)
  extra   = sorted(tips - taxa_set)
  print(f"[SKIP] {oma}: taxa mismatch TREE vs PHYLIP", file=sys.stderr)
  print(f"       PHYLIP({len(taxa_set)}): {','.join(sorted(taxa_set))}", file=sys.stderr)
  print(f"       TREE ({len(tips)}): {','.join(sorted(tips))}", file=sys.stderr)
  if missing: print(f"       missing_in_tree: {','.join(missing)}", file=sys.stderr)
  if extra:   print(f"       extra_in_tree: {','.join(extra)}", file=sys.stderr)
  sys.exit(25)

(outdir / f"{oma}.tree").write_text(tree + "\n")
PY
  then
    echo "⚠️  SKIP $oma: prep failed (see stderr above)"
    ((n_skip+=1))
    continue
  fi

  # codeml control
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

        icode = 4

   fix_kappa = 0
       kappa = 2.0

   fix_omega = 0
       omega = 1

  cleandata = 0
 RateAncestor = 2
EOF

  echo "▶ Running codeml for $oma ..."
  if ( cd "$outdir" && printf '\n' | codeml codeml.ctl > codeml.log 2>&1 ); then
    echo "✅ $oma done"
    ((n_ok+=1))
  else
    echo "❌ $oma FAILED (see $outdir/codeml.log)"
    ((n_fail+=1))
  fi
done

echo ""
echo "================ SUMMARY ================"
echo "OK:    $n_ok"
echo "FAIL:  $n_fail"
echo "SKIP:  $n_skip"
echo "TOTAL: ${#ALNS[@]}"
echo "========================================="

