#!/usr/bin/env bash
# 
# Plot IBD between each sample and its putative ancestors in windows across the
# genome.
#
# This runs `methlab ibdpainting` in a loop on each sample inside VCF files for
# each sequencing plate. Job arrays indicate which plate to run.
# Input:
#     # VCF files of genotype calls at SNPs segregating between the four parental
#         accessions.
# Output:
    # Output of ibdpainting for each sample, sorted by plate:
        # Plot of IBD across the genome
        # Scores for each candidate across windows.

# SLURM
#SBATCH --job-name=ibdpainting
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=2GB
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=1058-1153

# Load the conda environment
source setup.sh

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

indir=03_processing/04_validate_samples/output

# Sample sheet that includes ancestry
sample_sheet=$indir/sample_sheet.csv
# Information about sequencing plate, plantID, and ancestry.
plate=$(  awk -v row="$i" -F',' 'NR==row {print $1}' $sample_sheet)
plantID=$(awk -v row="$i" -F',' 'NR==row {print $2}' $sample_sheet)
crossID=$(awk -v row="$i" -F',' 'NR==row {print $6}' $sample_sheet)

# HDF5 file for the parents
parents_hdf5=$indir/parents_genic_SNPs_only.hdf5
# Directory containing the HDF5 files for the test individuals
test_dir=$indir/hdf5



# === Output files === #

outdir=$indir/ibdpainting/$plate
mkdir -p $outdir



# === Script === #

# If the sample is a Col0 positive control, set the crossID to look for this genome
if [[ $plantID == *Col0* ]]; then
    crossID = "6909"
fi

# Run the validation program.
ibdpainting \
    --input $test_dir/${plate}.hdf5 \
    --reference $parents_hdf5 \
    --sample_name $plantID \
    --window_size 500000 \
    --outdir $outdir \
    --expected_match ${crossID/x/ }