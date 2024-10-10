#!/usr/bin/env bash
# 
# For plate 2021-007 we repeated sequencing of the plate to get additional 
# coverage. This script merges fastq files from the two runs prior to alignment.
#
# Input:
    # Directories of unaligned fastq files.
# Output:
    # Directory of merged fastq files.

# SLURM
#SBATCH --job-name=merge_fastq
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=rapid
#SBATCH --time=5:00
#SBATCH --array=0-191


# Define working directory and load the conda environment
source setup.sh
i=${SLURM_ARRAY_TASK_ID}
# === Input paths ====

# Path to fastq files for the initial run 
indir=($workdir/04_validate_samples/01_raw_bs_reads/22H7YVLT3_3_R16762_20240222/demultiplexed/282462/*_R{1,2}_*)
# Path to the first file of the pair.
file1=${indir[$i]}
echo "Path to the first file of the pair:\n ${file1}\n"

# === Output paths ====

outdir=$workdir/04_validate_samples/02_merge_fastq
mkdir -p $outdir

# === Script ===

# Swap out the run information for the file in the first run
file2=${file1/22H7YVLT3_3_R16762_20240222/22KCNVLT3_7_R17286_20240606}
file2=${file2//282462/295386}
echo "Path to the second file of the pair:\n ${file2}\n"
# Swap request number for "merged" to create a name for the output file.
outfile=$(basename $file1)
outfile=${outfile/282462/merged}
echo "Saving merged fastq file to:\n ${outdir}/${outfile}\n"
# Merge fastq files.
cat $file1 $file2 > $outdir/$outfile