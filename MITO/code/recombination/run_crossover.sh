#!/bin/bash
#SBATCH --partition=general

#SBATCH --qos=grp_mlynch11

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --time=2-00:00:00

set -euo pipefail



# ============================================================

# LDhat pairwise — CROSSOVER ONLY (L.locs) — ALL SPECIES

# Runs BOTH theta values per species (theta_low + theta_high)

# Saves outputs under:

#   sites_locs/<species>/crossover/theta_<THETA>/

#     pairwise.stdout, pairwise.stderr, summary.txt

#

# Matches your exact interactive prompt order, INCLUDING:

#   "Do you wish to test for recombination? (yes=1, no=0): 1"

#

# Prompts answered in this order (from your run):

#   1) theta per site

#   2) Max 4Ner for grid (suggest 100)           -> 100

#   3) Number of points on grid (suggest 101)    -> 101

#   4) change grid? (1/0)                        -> 0

#   5) sliding windows? (first time)             -> 0

#   6) sliding windows? (second time)            -> 0

#   7) Rmin table? (0/1/2)                       -> 1

#   8) moment method? (1/0)                      -> 1

#   9) test for recombination? (1/0)             -> 1

# (pairwise then runs 1000 shuffles by default; no extra input needed)

# ============================================================



LDHAT_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/LDhat_install/LDhat"

PAIRWISE="${LDHAT_DIR}/pairwise"



SITESLOCS_ROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/sites_locs"

OUTROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/sites_locs"



RHOMAX=100

NPTS=101



# ---------- run list (species theta_low theta_high) ----------

run_species() {

  cat <<'EOF'

biaurelia        0.01    0.1

bursaria         0.001   0.01

caudatum         0.01    0.1

pfoki           0.01    0.1

multi            0.001   0.01

pquin            0.001   0.01

primaurelia      0.001   0.01

quadecaurelia    0.0001  0.001

sexaurelia       0.01    0.1

tetraurelia      0.001   0.01

triaurelia       0.001   0.01

tribursaria      0.001   0.01

EOF

}



# ---------- main ----------

while read -r sp th_low th_high; do

  [[ -z "${sp:-}" || "$sp" == \#* ]] && continue



  SPDIR="${SITESLOCS_ROOT}/${sp}"

  SITES="${SPDIR}/${sp}.sites"

  LOCL="${SPDIR}/${sp}.L.locs"



  [[ -f "$SITES" ]] || { echo "ERROR: missing $SITES"; exit 1; }

  [[ -f "$LOCL"  ]] || { echo "ERROR: missing $LOCL";  exit 1; }



  for THETA in "$th_low" "$th_high"; do

    OUTDIR="${OUTROOT}/${sp}/crossover/theta_${THETA}"

    mkdir -p "$OUTDIR"

    cd "$OUTDIR"



    echo "▶ ${sp} crossover theta=${THETA}"



    "${PAIRWISE}" -seq "$SITES" -loc "$LOCL" > pairwise.stdout 2> pairwise.stderr <<EOF

${THETA}

${RHOMAX}

${NPTS}

0

0

0

1

1

1

EOF



    # Extract the exact lines you care about (table-ready)

    awk '

      /Observed LD statistics/ ||

      /Npairs/ || /Lkmax/ || /G4/ || /corr\(r2/ || /corr\(D/ ||

      /Results of shuffling/ ||

      /Proportion Lkmax/ || /Proportion G4/ || /Proportion corr\(r2/ || /Proportion corr\(D/

    ' pairwise.stdout > summary.txt



    echo "✅ wrote: $OUTDIR/summary.txt"

  done

done < <(run_species)



echo "🎉 All crossover runs finished."

echo "See outputs under: ${OUTROOT}/<species>/crossover/theta_<THETA>/"


