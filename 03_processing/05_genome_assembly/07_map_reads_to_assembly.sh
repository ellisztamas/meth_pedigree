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
#     Fastq file containing raw reads
#     FASTA file for a scaffolded assembly to map to 
# Output:
#       
#
# Tom Ellis, 18th October 2024, adapting code by Robin Burns.

# SLURM
#SBATCH --job-name=maps_reads_to_assembly
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4:00:00
#SBATCH --array=0-3

# set -e

# Define working directory and load the conda environment
source setup.sh

# === Input === #

i=$SLURM_ARRAY_TASK_ID

# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Genome to map to
# scaffolded_assembly=$workdir/05_genome_assembly/03_scaffolding/${accession_name}/8334/${accession_name}_8334_scaffolded_contigs.fasta
ref_genome=01_data/04_parental_genomes/8334.fa

# Input raw reads
input_raw_reads=$workdir/05_genome_assembly/01_raw_data/${accession_name}.fastq.gz


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

# Text file giving contig, bp position and read depth
mosdepth_prefix=${outdir}/${accession_name}


# === Script === #

# Create a k-mer library for mapping
echo -n "Creating k-mer library..."
SECONDS=0
meryl count k=15 output $meryldb ${ref_genome}
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
    ${ref_genome} \
    ${input_raw_reads} > \
    $sam_alignment
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds.\n"

# Convert the SAM file to BAM and sort it 
echo -n "Sorting out the output bam alignment..."
SECONDS=0
samtools faidx $ref_genome
samtools view \
    -@ ${SLURM_CPUS_PER_TASK} -bh \
    -t ${ref_genome}.fai \
    -o $unsorted_bam \
    $sam_alignment
samtools sort \
    -@ ${SLURM_CPUS_PER_TASK} \
    -o ${sorted_bam} \
    ${unsorted_bam}
samtools index ${sorted_bam}
duration=$SECONDS
echo "finished in $((duration / 60)) minutes and $((duration % 60)) seconds.\n"

# Get coverage in windows of 1000 bp
mosdepth \
    --by 1000 \
    --no-per-base \
    -t ${SLURM_CPUS_PER_TASK} \
    $mosdepth_prefix \
    $sorted_bam