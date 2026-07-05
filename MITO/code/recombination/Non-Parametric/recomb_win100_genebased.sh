#!/bin/bash

#SBATCH -p general

#SBATCH -q grp_mlynch11

#SBATCH -J phipack_gene_based

#SBATCH -t 07:00:00

#SBATCH --mem=32G

#SBATCH --cpus-per-task=4

#SBATCH -o phipack_gene_based_%j.out

#SBATCH -e phipack_gene_based_%j.err



set -euo pipefail



module purge

module load mamba/latest

source activate pyenv



# ---- Paths ----

BASE_PN="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pn_ps_new"

OUTROOT="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/gene_based"



# Species directories present in pn_ps_new

SPECIES_DIRS=(

  biaurelia

  bursaria

  caudatum

  fokini

  primaurelia

  quadecaurelia

  sexaurelia

  tetraurelia

  tredecaurelia

  triaurelia

  tribursaria

  pquin

  multi

)



WINDOW=100



echo "[INFO] Starting PhiPack runs with window size ${WINDOW}"

echo "[INFO] Base pn_ps_new dir: ${BASE_PN}"

echo "[INFO] Output root: ${OUTROOT}"



for sp in "${SPECIES_DIRS[@]}"; do

    INDIR="${BASE_PN}/${sp}/codon_alignment"



    if [[ ! -d "$INDIR" ]]; then

        echo "⚠️  Skipping species '${sp}': directory not found: ${INDIR}"

        continue

    fi



    echo "====== Species: ${sp} ======"

    mkdir -p "${OUTROOT}/${sp}"



    shopt -s nullglob

    for f in "${INDIR}"/*_codon_aligned.fna; do

        [[ -e "$f" ]] || continue



        fname=$(basename "$f")

        # Strip the suffix "_codon_aligned.fna" to get the gene name

        gene=${fname%_codon_aligned.fna}



        outfile="${OUTROOT}/${sp}/${gene}_phi_w${WINDOW}.log"



        echo "🔍 Running PhiPack on ${f}"

        echo "   → ${outfile}"



        # --- basic QC: how many sequences? ---

        nseq=$(grep -c '^>' "$f" || true)



        if (( nseq < 2 )); then

            echo "Gene ${gene} in ${sp}: fewer than 2 sequences in alignment. Skipping Phi test and marking as absent." > "${outfile}"

            continue

        fi



        # Run PhiPack with window=100 and 1000 permutations.

        # If Phi fails (e.g., too few informative sites), we catch it and mark as absent.

        if ! Phi -f "$f" -w "${WINDOW}" -o -p 1000 > "${outfile}" 2>&1; then

            echo "" >> "${outfile}"

            echo "Phi failed for ${gene} in ${sp} (likely too few informative sites / too short alignment)." >> "${outfile}"

            echo "Treating Phi result as ABSENT for downstream analysis." >> "${outfile}"

        fi

    done

    shopt -u nullglob



    echo "✅ Completed species: ${sp}"

done



echo "🎉 All PhiPack runs finished (window=${WINDOW})."


