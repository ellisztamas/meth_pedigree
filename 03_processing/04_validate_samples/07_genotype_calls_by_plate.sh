#!/usr/bin/env bash

# Call SNPs from bisulphite data at sites polymorphic in the four parental
# accessions and within genes.
#
# For each plate separately, this uses `bcftools mpileup` to call SNPs at genic
# positions only. It tidies up the resulting VCF files, and converts them to 
# HDF5 files using the function `vcf_to_hdf5` in scikit-allel.
#
# Input:
#    aligned, sorted, depulicated BAM files
#    Reference FASTA genome
#    Text file giving SNP positions to call, in genes only. See 06_create_parental_VCF.sh
# Output:
#    HDF5 file of a VCF file for each plate.

# Tom Ellis, 19th June 2024

# SLURM
#SBATCH --job-name=genotype_calls
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=long
#SBATCH --time=4-00:00:00
#SBATCH --array=11

# Set working directory and load conda environment
source setup.sh
# set -e

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory containing aligned BAM files
indirs=($workdir/04_validate_samples/04_align_to_tair10/202?-0??)
plate=${indirs[$i]}
echo "Input plate = ${plate}."

# Location of the reference genome to map to.
genome=03_processing/01_prepare_genomes/tair10/TAIR10_plus_vectors.fa

# Table listing SNP positions that vary between the parents and are found in genes.
variable_SNP_positions=03_processing/04_validate_samples/output/variable_SNP_positions.tsv.gz

# === Output files === #
# Output directory
outdir=$workdir/04_validate_samples/07_genotype_calls_by_plate
mkdir -p $outdir

# file with location of bamfiles (one per line)
bam_list=${outdir}/bam_list_"$(basename ${plate})".txt

# Name for the output VCF file
vcf_with_full_paths=$outdir/$(basename ${plate})_full_paths.vcf.gz
vcf_with_sample_names=${vcf_with_full_paths/_full_paths/} # VCF with sample names renamed to exclude full paths.

# Files giving sample names before and after pruning the path name
old_names=${vcf_with_full_paths/full_paths.vcf.gz/old_names.txt}
new_names=${old_names/old_/new_}

# HDF5 file
proj_dir=03_processing/04_validate_samples/output/hdf5
mkdir -p $proj_dir
hdf5_file=$proj_dir/$(basename ${plate}).hdf5

# === Script === #

# create file with bamfile paths
ls -d ${plate}/aligned_bams/*.bam > $bam_list

# Run the SNP calls
echo "Getting genotype likelihoods, and using them to call SNPs"
if [ -f ${vcf_with_full_paths} ]; then rm $vcf_with_full_paths ; fi
bcftools mpileup --min-MQ 20 -a FORMAT/DP --skip-indels -f $genome -b $bam_list -Ou | \
bcftools call -m --constrain alleles --targets-file $variable_SNP_positions -Oz --output $vcf_with_full_paths

# Rename the samples to remove absolute paths and .bam suffixes
# Create text files giving the old and new names
echo "Renaming samples in the VCF file."
bcftools query -l $vcf_with_full_paths > $old_names
xargs -rd '\n' -a $old_names basename -a --suffix=.sortedByPos.bam > $new_names
# Swap the header_names
if [ -f $vcf_with_sample_names ]; then rm $vcf_with_sample_names; fi
bcftools reheader \
    -s $new_names \
    -o $vcf_with_sample_names \
    $vcf_with_full_paths
if [ $? -eq 0 ] ; then echo "Renamed samples successfully"; fi
# Index the VCF file
tabix $vcf_with_sample_names

# Convert VCF to HDF5.
echo "Creating HDF5 file."
02_library/vcf_to_HDF5.py --input $vcf_with_sample_names --output $hdf5_file
