#!/usr/bin/env bash
# 
# Assess the completeness of the assembly with Merqury
#
# Input:
#     FASTA file containing contigs that could be scaffolded to an autosome in
#          the reference genome.
# Output:
#    Meryl database files
#    Output of Merqury
#       
#
# Tom Ellis, 22nd November 2024

# SLURM
#SBATCH --job-name=merqury
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --array=0-3

# set -e

# Define working directory and load the conda environment
source setup.sh

# === Input === #

i=$SLURM_ARRAY_TASK_ID

# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Input raw fastq file
input_raw_reads=$workdir/05_genome_assembly/01_raw_data/${accession_name}.fastq.gz
echo "Processing reads from input file: ${input_raw_reads}"

# FASTA file containing contigs for the autosomes
scaffolded_assembly=${workdir}/05_genome_assembly/03_scaffolding/${accession_name}/8334/${accession_name}_scaffolded_contigs.fasta


# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/05_merqury/$accession_name
mkdir -p $outdir

# K-mer library files
meryldb=$outdir/${accession_name}.meryl

# Prefix for the Merqury output
output_prefix=$accession_name

# === Script === #

cd $outdir

# Create a k-mer library for mapping
echo -n "Creating k-mer library..."
SECONDS=0
meryl k=18 count $scaffolded_assembly output $meryldb
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds."

# Run Merqury
echo -n "Analysing the genome with Merqury..."
SECONDS=0
$MERQURY/merqury.sh \
    $meryldb \
    $scaffolded_assembly \
    $output_prefix
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds."
