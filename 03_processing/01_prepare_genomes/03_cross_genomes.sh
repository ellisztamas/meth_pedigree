#!/bin/bash

# Script to merge genome assemblies for the parents of the crosses with the one another,
# the TAIR10 choroplast and mitochondria, as well as the lambda, phi X and pUC19 controls.
# The assembly for 6024 already has the organelles, so these are not added again.
#
# The merged fasta files are then prepared for bisulphite read mapping by running 
# them bismark_genome_preparation.
#
# Input:
#     Fasta files for parents, organelles, lambda phage and pUC19 vectors
# Output: 
#     A concatenated fasta file for each of four crosses in its own directory.
#     Bisulphite-converted genome files for each cross to use with Bismark
#
# Tom Ellis, 25th March 2024

#SBATCH --job-name=prepare_cross_genomes
#SBATCH --time=3:00:00
#SBATCH --qos=short
#SBATCH --mem=10gb
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err

date

# Load the conda environment
source setup.sh

# source directories containing parental and vector genomes
parents_dir=01_data/04_parental_genomes
tair10_dir=01_data/03_tair10
vector_dir=01_data/07_vector_DNA

# 1158x6184
outdir=03_processing/01_prepare_genomes/1158x6184
mkdir -p $outdir
# Concatenate the autosome assembly, organelles and vectors
cat $parents_dir/1158_v1.1.fasta \
    $parents_dir/6184_v1.1.fasta \
    $vector_dir/pUC19.fasta \
    $vector_dir/lambda_vector.fasta \
    $vector_dir/phiX.fasta \
    $tair10_dir/TAIR10_chrC.fas \
    $tair10_dir/TAIR10_chrM.fas > \
    $outdir/1158x6184.fa
# Bisulphite convert the genome
bismark_genome_preparation $outdir
if [ $? -eq 0 ] ; then echo "Successfully created bisulphite-converted genome for 1158x6184" ; fi

# 1158x6024
outdir=03_processing/01_prepare_genomes/1158x6024
mkdir -p $outdir
# Concatenate the assembly and vectors
cat $parents_dir/1158_v1.1.fasta \
    $parents_dir/AT_6024.v2.1.addChrMC.fasta \
    $vector_dir/pUC19.fasta \
    $vector_dir/lambda_vector.fasta \
    $vector_dir/phiX.fasta > \
    $outdir/1158x6024.fa
# Bisulphite convert the genome
bismark_genome_preparation $outdir
if [ $? -eq 0 ] ; then echo "Successfully created bisulphite-converted genome for 1158x6024" ; fi

# 1158x8249
outdir=03_processing/01_prepare_genomes/1158x8249
mkdir -p $outdir
# Concatenate the autosome assembly, organelles and vectors
cat $parents_dir/1158_v1.1.fasta \
    $parents_dir/8249_v1.1.fasta \
    $vector_dir/pUC19.fasta \
    $vector_dir/lambda_vector.fasta \
    $vector_dir/phiX.fasta \
    $tair10_dir/TAIR10_chrC.fas \
    $tair10_dir/TAIR10_chrM.fas > \
    $outdir/1158x8249.fa
# Bisulphite convert the genome
bismark_genome_preparation $outdir
if [ $? -eq 0 ] ; then echo "Successfully created bisulphite-converted genome for 1158x8249" ; fi


# 8249x6024
outdir=03_processing/01_prepare_genomes/8249x6024
mkdir -p $outdir
# Concatenate the assembly and vectors
cat $parents_dir/8249_v1.1.fasta \
    $parents_dir/AT_6024.v2.1.addChrMC.fasta \
    $vector_dir/pUC19.fasta \
    $vector_dir/lambda_vector.fasta \
    $vector_dir/phiX.fasta > \
    $outdir/8249x6024.fa
# Bisulphite convert the genome
bismark_genome_preparation $outdir
if [ $? -eq 0 ] ; then echo "Successfully created bisulphite-converted genome for 8249x6024" ; fi



# 8249x6184
outdir=03_processing/01_prepare_genomes/8249x6184
mkdir -p $outdir
# Concatenate the autosome assembly, organelles and vectors
cat $parents_dir/8249_v1.1.fasta \
    $parents_dir/6184_v1.1.fasta \
    $vector_dir/pUC19.fasta \
    $vector_dir/lambda_vector.fasta \
    $vector_dir/phiX.fasta \
    $tair10_dir/TAIR10_chrC.fas \
    $tair10_dir/TAIR10_chrM.fas > \
    $outdir/8249x6184.fa
# Bisulphite convert the genome
bismark_genome_preparation $outdir
if [ $? -eq 0 ] ; then echo "Successfully created bisulphite-converted genome for 8249x6184" ; fi

# 6184x6024
outdir=03_processing/01_prepare_genomes/6184x6024
mkdir -p $outdir
# Concatenate the assembly and vectors
cat $parents_dir/6184_v1.1.fasta \
    $parents_dir/AT_6024.v2.1.addChrMC.fasta \
    $vector_dir/pUC19.fasta \
    $vector_dir/lambda_vector.fasta \
    $vector_dir/phiX.fasta > \
    $outdir/6184x6024.fa
# Bisulphite convert the genome
bismark_genome_preparation $outdir
if [ $? -eq 0 ] ; then echo "Successfully created bisulphite-converted genome for 6184x6024" ; fi

date