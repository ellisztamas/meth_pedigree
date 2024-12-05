#!/usr/bin/env bash
# 
# Script to convert raw BAM files to fastq, then assemble with hifiasm.
# Unscaffolded contigs are converted to FASTA.
#
# Input:
#       Raw FASTQ files for each sample
# Output:
#       Hifi assembly files. We don't need phased assemblies, so the main
#            assembly file is the most interesting (ends in .bp.p_ctg.gfa)
#       FASTA file containing unscaffolded contigs
#
# Tom Ellis, 18th October 2024, adapting code by Robin Burns.

# SLURM
#SBATCH --job-name=hifi_assembly
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=38
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --array=2-3

set -e

# Define working directory and load the conda environment
source setup.sh

# === Input === #

i=$SLURM_ARRAY_TASK_ID

# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Input raw fastq file
raw_read_array=($workdir/05_genome_assembly/01_raw_data/${accession_name}.fastq.gz)
input_raw_reads=${raw_read_array[$i]}
echo "Processing reads from input file: ${raw_bam}"


# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/02_assembly
mkdir -p $outdir

# Prefix for various output files is the directory and the accession name
hifiasm_prefix=$outdir/$accession_name

# Contig FASTA file
unscaffolded_contigs=${outdir}/${accession_name}_unscaffolded_contigs.fa



# === Script === #

# Assembly
echo -n "Beginning the assembly with hifiasm..."
SECONDS=0
hifiasm \
    -o $hifiasm_prefix \
    -t ${SLURM_CPUS_PER_TASK} \
    -f 0 \
    -l 0 \
    ${input_raw_reads}
duration=$SECONDS
echo "finished in $((duration / 3600)) hours and $((duration / 60)) minutes.\n"

# Create a FASTA file for the unscaffolded contigs
echo "Converting GFA to FASTA"
awk '/^S/{print ">"$2; print $3}' ${hifiasm_prefix}.bp.p_ctg.gfa > ${unscaffolded_contigs}