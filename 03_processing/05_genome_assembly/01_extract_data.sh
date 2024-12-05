#!/usr/bin/env bash
# 
# Unpack the data to scratch-cbe
#
# Input:
#       Tarball containing raw reads
# Output:
#       Unzipped tarball
#
# Tom Ellis, 21st October 2024.

# SLURM
#SBATCH --job-name=extract_data
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --qos=short
#SBATCH --mem=5G
#SBATCH --time=8:00:00

# Define working directory and load the conda environment
source setup.sh

# === Input === #

raw_reads=01_data/09_pacbio_reads/slskuGak5s-data.tar
# raw_reads=/groups/nordborg/raw.data/athaliana/dna/PacBio/R17196/slskuGak5s-data.tar

# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/01_raw_data
# outdir=01_data/09_pacbio_reads/data/hifi_reads
mkdir -p $outdir

# === Script === #

# Unpacking the raw data
tar -xvC $outdir -f $raw_reads

# Convert bam to fastq
bam2fastq -o $outdir/1158 $outdir/data/hifi_reads/m84151_240831_053551_s2.hifi_reads.bc2077.bam
bam2fastq -o $outdir/6184 $outdir/data/hifi_reads/m84151_240831_053551_s2.hifi_reads.bc2078.bam
bam2fastq -o $outdir/8249 $outdir/data/hifi_reads/m84151_240831_053551_s2.hifi_reads.bc2079.bam
bam2fastq -o $outdir/6024 $outdir/data/hifi_reads/m84151_240831_053551_s2.hifi_reads.bc2080.bam

