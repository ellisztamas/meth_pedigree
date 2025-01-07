#!/usr/bin/env bash
# 
# Concatenate assembled autosomes with TAIR10 organelles.
#
# Tom Ellis, 7th January 2025

# SLURM
#SBATCH --job-name=cat_organelles
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=rapid
#SBATCH --time=10:00
#SBATCH --array=0-3

# Define working directory and load the conda environment
source setup.sh

# === Input === #

i=${SLURM_ARRAY_TASK_ID}

# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Autosome assemblies
# Use the main assemblies scaffolded to 8334 for three accessions
# For 6024, use the assembly using the split contig ptg000002l
asm_array=(
    03_processing/05_genome_assembly/output/03_scaffolding/1158/8334/1158_8334_scaffolded_contigs.fasta
    03_processing/05_genome_assembly/output/03_scaffolding/6184/8334/6184_8334_scaffolded_contigs.fasta
    03_processing/05_genome_assembly/output/03_scaffolding/8249/8334/8249_8334_scaffolded_contigs.fasta
    03_processing/05_genome_assembly/output/11_fix_6024/6024_8334_scaffolded_contigs.fasta
    )

# Organelles
# I will also define these as arrays even though they are all from TAIR10 in case
# I manage to get individual assemblies later
chloroplast_array=(
    01_data/03_tair10/TAIR10_chrC.fas
    01_data/03_tair10/TAIR10_chrC.fas
    01_data/03_tair10/TAIR10_chrC.fas
    01_data/03_tair10/TAIR10_chrC.fas
)
mitochondrion_array=(
    01_data/03_tair10/TAIR10_chrM.fas
    01_data/03_tair10/TAIR10_chrM.fas
    01_data/03_tair10/TAIR10_chrM.fas
    01_data/03_tair10/TAIR10_chrM.fas
)



# === Output === #

# Directory to store the output
outdir=${workdir}/05_genome_assembly/12_concatenate_organelles
mkdir -p $outdir

# Autosomes with headers Chr1, Chr2 etc (removing "_RagTag")
corrected_headers=$outdir/${accession_name}_autosomes.fasta
# Autosomes plus organelles
outfile=$outdir/${accession_name}.fasta

# Project directory
proj_dir=03_processing/05_genome_assembly/output/


# === Script === #

# Remove "_RagTag" suffixes from FASTA headers
awk '/^>/ {sub("_RagTag", "", $1)} 1' ${asm_array[$i]} > $corrected_headers

# Concatenate autosomes and organelles
cat $corrected_headers ${chloroplast_array[$i]} ${mitochondrion_array[$i]} > $outfile

# Stage out
cp $outfile $proj_dir