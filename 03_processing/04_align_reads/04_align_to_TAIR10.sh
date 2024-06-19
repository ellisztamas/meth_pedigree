#!/bin/bash

# Align raw reads to the TAIR10 reference genome.
#
# Note that the script hard codes various optional arguments to trim galore and
# Bismark - check that script for which arguments are used.
#
# Input:
#     Directories containing zipped fastq files
#     CSV sample sheet giving
#        - plantID
#        - paths to raw fastq files
#        - paths to genome to assemble to 
# Output:
#   - logs
#     - bismark_report.txt   
#     - bismark_deduplication_report.txt
#   - reports
#     - CX reports for each samples
#   - sorted
#       Aligned BAM files sorted by position
#       Index files for those BAM files
#
# Tom Ellis, 20th March 2024

#SBATCH --job-name=align_tair10
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --time=16:00:00
#SBATCH --qos=medium
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5G
#SBATCH --array=960-961

# Set up working directory and load conda evironment
source setup.sh

i=${SLURM_ARRAY_TASK_ID}

# === Input files === #
# Sample sheet giving sample name and paths to the two fastq files
sample_sheet=03_processing/04_align_reads/output/sample_sheet.csv
# Path to the assembly to map to
genome=03_processing/01_prepare_genomes/tair10/TAIR10_plus_vectors.fa
# Input settings to run the script
trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 9 --three_prime_clip_R2 9"
bismark_args="--local --non_directional --strandID"
methylation_extractor_args="--cytosine_report --CX_context --no_header --no_overlap --comprehensive --gzip"

# === Output files === #
# Directory to save the main output of the pipeline
scratch_dir=$workdir/04_align_reads/03_align_to_tair10
mkdir -p $scratch_dir
# Directory to copy the really interesting files to keep
proj_dir=03_processing/04_align_reads/output/tair10

# === Script === #

# Get the plate ID
plate_col=$(cut -d',' -f 1 $sample_sheet)
plate_col=($plate_col)
plate=${plate_col[$i]}
mkdir -p $proj_dir/${plate} # Create the directory for the output
# Get the sample name
sample_names=$(cut -d',' -f 2 $sample_sheet)
sample_names=($sample_names)
plantID=${sample_names[$i]}
# Path for read 1
read1_col=$(cut -d',' -f 3 $sample_sheet)
read1_col=($read1_col)
read1=${read1_col[$i]}
# Path for read 2
read2_col=$(cut -d',' -f 4 $sample_sheet)
read2_col=($read2_col)
read2=${read2_col[$i]}

# Run the pipeline
02_library/bismark_pipeline.sh \
    --sample $plantID \
    --read1  $read1 \
    --read2  $read2 \
    --genome $genome \
    --work $scratch_dir \
    --outdir $proj_dir/${plate} \
    --trim_galore_args "${trim_galore_args}" \
    --bismark_args "${bismark_args}" \
    --methylation_extractor_args "${methylation_extractor_args}"



# if [ $? -eq 0 ] ; then echo -n "Bismark script run successfully" ; fi

date