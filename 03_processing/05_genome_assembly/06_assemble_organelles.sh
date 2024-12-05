#!/usr/bin/env bash
# 
# Assemble chloroplast and mitochondria using TIPP.
#
# Note that this uses a different conda environment from 
# Input:
#     Fastq files containing raw PacBio reads
# Output:
#     
#
# Tom Ellis, 25th November 2024

# SLURM
#SBATCH --job-name=organelle_assembly
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --array=0

set -e

# Define working directory and load the conda environment
# (because this is different from what's in setup.sh)
workdir=/scratch-cbe/users/$(whoami)/meth_pedigree
# source activate /users/$(whoami)/.conda/envs/TIPP

# === Input === #

i=$SLURM_ARRAY_TASK_ID
# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Input raw fastq file
raw_read_array=($workdir/05_genome_assembly/01_raw_data/${accession_name}.fastq.gz)
input_raw_reads=${raw_read_array[$i]}


# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/06_assemble_organelles/$accession_name
mkdir -p $outdir


# === Script === #

cd $outdir

apptainer pull docker://weigelworld/tipp:latest

# TIPP_plastid.v2.1.pl -f $input_raw_reads -t ${SLURM_CPUS_PER_TASK}
singularity exec tipp_latest.sif /usr/local/bin/TIPP_plastid.v2.1.pl -f $input_raw_reads -t ${SLURM_CPUS_PER_TASK}