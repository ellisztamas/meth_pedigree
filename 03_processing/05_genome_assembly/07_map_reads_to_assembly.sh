#!/usr/bin/env bash
# 
# Map raw pacbio reads to the pacbio assembly.
#
# This uses meryl to count kmers, and feeds those into winnowmap to align reads
# to the genome.
# 
# At present this uses the assembly scaffolded to TAIR10. This needs updating
# to use an assembly from 03_scaffolding. 
#
# Input:
#       Raw BAM files for each sample
# Output:
#       Zipped fastq file
#       Hifi assembly files. We don't need phased assemblies, so the main
#            assembly file is the most interesting (ends in .bp.p_ctg.gfa)
#       FASTA file for the assembly
#       Output of ragtag.
#
# Tom Ellis, 18th October 2024, adapting code by Robin Burns.

# SLURM
#SBATCH --job-name=maps_reads_to_assembly
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G
#SBATCH --time=24:00:00
#SBATCH --array=0-3

# set -e

# Define working directory and load the conda environment
source setup.sh

# === Input === #

i=$SLURM_ARRAY_TASK_ID

# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Input raw bam file
scaffolded_assembly=$workdir/05_genome_assembly/03_scaffolding/${accession_name}/8334/${accession_name}_8334_scaffolded_contigs.fasta

# Input raw reads
raw_read_array=($workdir/05_genome_assembly/01_raw_data/${accession_name}.fastq.gz)
input_raw_reads=${raw_read_array[$i]}

# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/07_map_reads_to_assembly/$accession_name
mkdir -p $outdir

# K-mer library files
meryldb=$outdir/meryldb
repetitive_k15=$outdir/repetitive_k15.txt

# Aligned SAM and BAM files
sam_alignment=${outdir}/${accession_name}.unsorted.sam
unsorted_bam=${sam_alignment/sam/bam}
sorted_bam=${unsorted_bam/unsorted/sorted}

# === Script === #

# Create a k-mer library for mapping
echo -n "Creating k-mer library..."
SECONDS=0
meryl count k=15 output $meryldb ${scaffolded_assembly}
meryl print greater-than distinct=0.9998 $meryldb > $repetitive_k15
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds.\n"

# Map raw reads to the assembly
echo -n "Mapping raw reads to the assembly..."
SECONDS=0
winnowmap \
    -W ${repetitive_k15} \
    -ax map-pb \
    -t ${SLURM_CPUS_PER_TASK} \
    ${scaffolded_assembly} \
    ${input_raw_reads} > \
    $sam_alignment
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds.\n"

# Convert the SAM file to BAM and sort it 
echo -n "Sorting out the output bam alignment..."
SECONDS=0
samtools faidx $scaffolded_assembly
samtools view \
    -@ ${SLURM_CPUS_PER_TASK} -bh \
    -t ${scaffolded_assembly}.fai \
    -o $unsorted_bam \
    $sam_alignment
samtools sort \
    -@ ${SLURM_CPUS_PER_TASK} \
    -o ${sorted_bam} \
    ${unsorted_bam}
samtools index ${sorted_bam}
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds.\n"
