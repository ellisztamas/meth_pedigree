#!/bin/bash

# Script to merge the TAIR10 fasta file with the lambda and pUC19 genomes
# This also creates a file with just the organelles for masking PacBio reads
#
# Input:
#     Fasta files for TAIR10, lambda phage and pUC19 vectors
# Output: 
#     A concatenated fasta file with all the inputs above
#     A concatenated fasta file with only the organelles
#     Bisulphite-converted genome files to use with Bismark
#
# Tom Ellis, 31st January 2024

#SBATCH --job-name=prepare_tair10_genome
#SBATCH --time=20:00
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err

date

# Load the conda environment
source setup.sh

outdir=03_processing/01_prepare_genomes/tair10
mkdir -p $outdir

# Concatenate the TAIR10 and Lambda phage fasta files
cat 01_data/03_tair10/TAIR10_chr_all.fas \
    01_data/07_vector_DNA/pUC19.fasta \
    01_data/07_vector_DNA/lambda_vector.fasta > \
    ${outdir}/TAIR10_plus_vectors.fa

# Concatenate the organelles only.
cat 01_data/03_tair10/TAIR10_chrC.fas \
    01_data/03_tair10/TAIR10_chrM.fas > \
    ${outdir}/tair10_organelles_only.fa


bismark_genome_preparation $outdir

date