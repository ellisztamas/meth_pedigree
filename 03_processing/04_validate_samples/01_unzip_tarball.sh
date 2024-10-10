#!/usr/bin/env bash
# 
# Script to unzip each raw data file to the `scratch-cbe` drive
# on the VBC cluster.
#
# Input:
    # Zipped tarball of fastq files
# Output:
    # Lots of unzipped files, including demultiplexed fastq folders
    # Directory of QC files copied to the project folder

# SLURM
#SBATCH --job-name=unzip_raw_reads
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=2GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-18

# Define working directory and load the conda environment
source setup.sh

# Directory with subdirectories for each plate
zipfile_array=(01_data/06_raw_bisulphite_reads/**/*.tar.gz)
# Location of the raw zip file on the VBC cluster
zipfile=${zipfile_array[$SLURM_ARRAY_TASK_ID]}

# Define a directory to stored upzipped files
scratchdir=$workdir/04_validate_samples/01_raw_bs_reads
# scratchdir=temp_unzip
mkdir -p $scratchdir
# Directory to copy the QC output
projdir=03_processing/04_validate_samples/output/multiqc_raw_reads
mkdir -p $projdir

# Unzip raw data
tar -xvC $scratchdir -f $zipfile

# copy the multiqc report to the project folder
cp -r $scratchdir/$(basename -s.tar.gz $zipfile)/qc/*_multiqc_report.html $projdir
