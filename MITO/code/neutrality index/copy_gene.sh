#!/usr/bin/env bash

set -euo pipefail



SRC="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/dn_ds/ancestral_dnds"

DEST="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/dn_ds/ancestral_dnds_filterallsites"



mkdir -p "$DEST"

shopt -s nullglob



for d in "$SRC"/*/; do

  gene="$(basename "$d")"

  f="$d/all_sequences_cleaned_with_gaps.fasta"

  if [[ -f "$f" ]]; then

    mkdir -p "$DEST/$gene"

    cp -vp "$f" "$DEST/$gene/${gene}_all_sequences_cleaned_with_gaps.fasta"

  else

    echo "Skipping $gene (no all_sequences_cleaned_with_gaps.fasta)" >&2

  fi

done


