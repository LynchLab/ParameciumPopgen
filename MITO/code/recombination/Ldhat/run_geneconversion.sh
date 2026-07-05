#!/bin/bash
#SBATCH --partition=general

#SBATCH --qos=grp_mlynch11

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --time=2-00:00:00

set -euo pipefail



# ============================================================

# LDhat pairwise — GENE CONVERSION ONLY (C.locs)

# Runs BOTH theta values per species (theta_low + theta_high)

# Saves outputs under:

#   <OUTROOT>/<species>/gene_conversion/theta_<THETA>/

#     pairwise.stdout, pairwise.stderr, summary.txt

#

# Uses the exact interactive answers you used:

#   tractlen=500

#   theta=<THETA>

#   rhomax=100

#   npts=101

#   change grid? 0

#   sliding windows? 0

#   Rmin table? 1  (Total only)

#   moment method? 1

#   test recombination? 1

# ============================================================



LDHAT_DIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/LDhat_install/LDhat"

PAIRWISE="${LDHAT_DIR}/pairwise"



SITESLOCS_ROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/sites_locs"

OUTROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/sites_locs"



TRACTLEN=500

RHOMAX=100

NPTS=101

NPERM=1000  # pairwise default shuffles in your run



# species  theta_low  theta_high

while read -r sp th_low th_high; do

  [[ -z "${sp:-}" || "$sp" == \#* ]] && continue



  SPDIR="${SITESLOCS_ROOT}/${sp}"

  SITES="${SPDIR}/${sp}.sites"

  LOCC="${SPDIR}/${sp}.C.locs"



  [[ -f "$SITES" ]] || { echo "ERROR: missing $SITES"; exit 1; }

  [[ -f "$LOCC"  ]] || { echo "ERROR: missing $LOCC";  exit 1; }



  for THETA in "$th_low" "$th_high"; do

    OUTDIR="${OUTROOT}/${sp}/gene_conversion/theta_${THETA}"

    mkdir -p "$OUTDIR"

    cd "$OUTDIR"



    echo "▶ ${sp} gene conversion theta=${THETA} tractlen=${TRACTLEN}"



    "${PAIRWISE}" -seq "$SITES" -loc "$LOCC" > pairwise.stdout 2> pairwise.stderr <<EOF

${TRACTLEN}

${THETA}

${RHOMAX}

${NPTS}

0

0

1

1

1

EOF



    # extract key stats + p-values (table-ready)

    awk '

      /Npairs/ || /Lkmax/ || /G4/ || /corr\(r2/ || /corr\(D/ ||

      /Proportion Lkmax/ || /Proportion G4/ || /Proportion corr\(r2/ || /Proportion corr\(D/

    ' pairwise.stdout > summary.txt



    echo "✅ wrote: $OUTDIR"

  done



done <<'EOF'

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



echo "🎉 All gene-conversion runs finished."

echo "See outputs under: ${OUTROOT}/<species>/gene_conversion/theta_<THETA>/"


