#!/bin/bash

#SBATCH --partition=general

#SBATCH --qos=grp_mlynch11

#SBATCH -J raxml_per_gene_copy_clean

#SBATCH -N 1

#SBATCH --mem=128G

#SBATCH --cpus-per-task=8

#SBATCH -t 7-00:00:00

#SBATCH -o raxml_per_gene_copy_clean_%j.out

#SBATCH -e raxml_per_gene_copy_clean_%j.err

#SBATCH --mail-user=$USER@asu.edu

#SBATCH --mail-type=BEGIN,END,FAIL



module purge

module load raxml-ng-1.1.0-openmpi-np



# Paths

BASEDIR="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/phylogeny_with_outgroup"

GENE_FASTA_DIR="${BASEDIR}/aligned"

OUTBASE="${BASEDIR}/gene_phylogeny"



mkdir -p "${OUTBASE}"



echo "Starting per-gene: COPY → CLEAN → RAxML..."

echo



for faa in "${GENE_FASTA_DIR}"/*.faa; do

    gene=$(basename "${faa}" .faa)



    echo "=== Processing gene: ${gene} ==="



    # Create gene output folder

    OUTDIR="${OUTBASE}/${gene}"

    mkdir -p "${OUTDIR}"



    # Step 1 — COPY original file into results folder

    RAW_COPY="${OUTDIR}/${gene}.faa"

    cp "${faa}" "${RAW_COPY}"

    echo "  Copied original FASTA → ${RAW_COPY}"



    # Step 2 — CLEAN the copied file

    CLEAN_MSA="${OUTDIR}/${gene}_clean.faa"

    echo "  Cleaning headers → ${CLEAN_MSA}"



    awk 'BEGIN{FS=OFS=""}

         /^>/{

            hdr=$0

            sub(/^>/,"",hdr)

            gsub(/[^A-Za-z0-9_.-]/,"_",hdr)

            print ">"hdr

            next

         }

         {print}' "${RAW_COPY}" > "${CLEAN_MSA}"



    # Show sample header before/after

    echo "  Sample header changes:"

    paste <(grep "^>" "${RAW_COPY}" | head -n 2) \

          <(grep "^>" "${CLEAN_MSA}" | head -n 2)

    echo



    # Step 3 — Run RAxML on the CLEANED FASTA

    PREFIX="${OUTDIR}/${gene}_raxml"



    raxml-ng-mpi --all --msa "${CLEAN_MSA}" --msa-format fasta --model LG+G --data-type AA --prefix "${PREFIX}" --threads 8 --tree rand{10} --bs-trees 10


    echo "  ✅ Finished gene: ${gene}"

    echo

done



echo "🎉 ALL DONE — Results in ${OUTBASE}/<gene>/"


