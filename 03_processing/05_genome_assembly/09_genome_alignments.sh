#!/usr/bin/env bash
# 
# Use pannagram to perform whole-genome alignment of the assemblies for 1158,
# 6024, 6184 an 8249 (scaffolded to TAIR12) and TAIR10. The reason for this is
# that scaffolding was better for TAIR12, but the annotations are better for
# TAIR10.
#
# Input:
#     FASTA file for a reference genome.
#     FASTA files for each genome assembly
# Output:
#     
#       
#
# Tom Ellis, 4th December 2024.

# SLURM
#SBATCH --job-name=ms_alignments
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --time=12:00:00

# Define working directory and load the conda environment
# (because this is different from what's in setup.sh)
workdir=/scratch-cbe/users/$(whoami)/meth_pedigree
source activate /users/$(whoami)/.conda/envs/pannagram


# === Input === #

# Reference genome
tair10=01_data/03_tair10/TAIR10_chr_all.fas

# Paths to scaffolded FASTA files.
indir=03_processing/05_genome_assembly/output/03_scaffolding

# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/09_genome_alignment

# Directory to copy the genomes
data_dir=$outdir/fasta_files
mkdir -p $data_dir


# === Script === #

# Pannagram needs the FASTA files to be in one directory.
# Copy the FASTA files to the working directory
cp $indir/**/tair10/*scaffolded_contigs.fasta $data_dir
cp $tair10 $data_dir

# Run the alignment
pannagram \
    -path_in $data_dir \
    -path_out $outdir \
    -cores ${SLURM_CPUS_PER_TASK} \
    -nchr 5 \
    -pre