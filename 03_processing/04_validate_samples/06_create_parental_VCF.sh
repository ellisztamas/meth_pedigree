#!/usr/bin/env bash

# Subset the published VCF from the 1001 genomes paper to include only the
# parental accessions and keep only variant sites within genes. Also create an
# Hdf5 version of the VCF file.
#
# Inputs:
#    VCF file from the 1001 genomes project
# Output:
#    Bed file giving start and end positions of annotaetd genes
#    Text file giving positions of SNPs to use to call SNPs in the crosses.
#    VCF file with four accessions, and only variable sites within genes.
#    HDF5 file from the VCF file

#
# Tom Ellis, 19th june 2024

# SLURM
#SBATCH --job-name=parental_vcf
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=5GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

# Set working directory
source setup.sh
set -e

# === Input files === #

# VCF file for 1163 inbred lines
input_VCF=01_data/08_1001genomes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz

# Annotated positions of genes in TAIR10.
genes_gff=01_data/03_tair10/TAIR10_GFF3_genes.gff

# === Output files === #

# Directory for output files.
outdir=$workdir/04_validate_samples/06_create_parental_vcf
mkdir -p $outdir

# BED file version of the GFF file
genes_bed=$outdir/TAIR10_GFF3_genes.bed

# VCF files for the four parental accessions
parents_VCF=$outdir/parents_only.vcf.gz

# File giving the positions of variable SNPs to use later
variable_SNP_positions=03_processing/04_validate_samples/output/variable_SNP_positions.tsv.gz

# HDF5 version of the parental VCF file
parents_hdf5=03_processing/04_validate_samples/output/parents_genic_SNPs_only.hdf5

# === Script === #

# Select only features labelled as coding genes in the GFF, and convert to BED
grep "Note=protein_coding_gene" $genes_gff | 
    awk -F'\t' '{print $1"\t"$4"\t"$5}' > $genes_bed

# Subset VCF for the four parental accessions, SNPs within genes that are present
# in at least on of the parents.
echo "Subsetting the parental VCF file."
bcftools view \
    --samples 1158,6024,6184,8249,6909 \
    -R $genes_bed \
    --min-ac 2 \
    --max-ac 8 \
    --output $parents_VCF \
    $input_VCF
tabix $parents_VCF

# Create a list of SNP positions to call in the experimental individuals.
echo "Creating targets file."
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $parents_VCF | \
    awk '{print $1, $2, $3}' OFS='\t' | \
    bgzip -c > $variable_SNP_positions
tabix -s1 -b2 -e2 $variable_SNP_positions

# Create an HDF5 file for the reference accessions
echo "Creating HDF5 file."
02_library/vcf_to_HDF5.py --input $parents_VCF --output $parents_hdf5
