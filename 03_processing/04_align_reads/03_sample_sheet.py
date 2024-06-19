# Create a sample sheet for aligning reads
#
# This looks up the plate, row and column position for a list of raw fastq files
# based on the barcodes in their filenames using methlab. These are merged with
# the real plantID for each individual. Finally this looks up the crossID from a
# file giving pedigree information.
#
# This script will throw a `SettingWithCopyWarning` for. As far as I can see the
# code is doing this correctly. This is just a warning and the code will
# probably run correctly. 
# 
# Input:
#     - CSV file giving the plantID and cross ID of every individual in the sample
#     - CSV file giving the plate, row and column position of each individual
#     - list of paths to raw fastq files.
# Output:
#     - CSV file with a row for each individual giving:
#         plate
#         plantID
#         paths to raw fastq files
#         path to the genome assembly to align reads to
#
# Tom Ellis, 18th March 2024

import pandas as pd
from glob import glob
import os

import methlab as ml
print("Using methlab version " + ml.__version__)

# CSV giving the cross ID of every individual in the sample
pedigree = pd.read_csv("03_processing/03_process_pedigree/output/pedigree_with_ancestry.csv")
# CSV file giving the plate, row and column position of each individual
ngs_sample_sheet = pd.read_csv("01_data/06_raw_bisulphite_reads/NGS_sample_sheet.csv", dtype=str)
# Directory containing raw reads.
scratchdir="/scratch-cbe/users/" + os.getlogin() + "/meth_pedigree/04_align_reads/01_raw_bs_reads/"
# scratchdir="/scratch-cbe/users/thomas.ellis/meth_pedigree/04_align_reads/01_raw_bs_reads/"
# scratchdir='temp_unzip/'

# ngs_fastq_sheet = ngs_sample_sheet.rename(columns={"sample":"plantID"})


# Dictionary of paths to the directory for each plate.
# 2021-007 is different because there were two sequencing runs, and the results
# FASTQ files merged into one.
fastqdir = {
    '2021-007' : '/scratch-cbe/users/thomas.ellis/meth_pedigree/04_align_reads/02_merge_fastq/',
    '2021-013' : scratchdir + '22KKY7LT3_2_R17070_20240503/demultiplexed/289773/',
    '2021-014' : scratchdir + '22KKY7LT3_3_R17071_20240430/demultiplexed/289774/',
    '2021-015' : scratchdir + '22KKCYLT3_4_R17029_20240413/demultiplexed/288527/',
    '2021-016' : scratchdir + '22KKCYLT3_3_R17028_20240413/demultiplexed/288526/',
    '2021-018' : scratchdir + '22KL2VLT3_2_R17129_20240509/demultiplexed/291136/',
    '2021-021' : scratchdir + '22KCNVLT3_5_R17269_20240606/demultiplexed/295267',
    '2021-023' : scratchdir + '22KL2VLT3_3_R17176_20240509/demultiplexed/293536/',
    '2022-005' : scratchdir + '22HC53LT3_5_R16842_20240306/demultiplexed/280348/',
    '2022-006' : scratchdir + '22HC53LT3_4_R16840_20240306/demultiplexed/280344/'
}

# Loop over each plate and line up file names with plate positions
fastq_sheet = {}
for k in fastqdir.keys():
    if os.path.exists(fastqdir[k]):
        mate1 = glob(fastqdir[k] + "/*_R1_*fastq.gz")
        mate2 = glob(fastqdir[k] + "/*_R2_*fastq.gz")
        sample_sheet_plate_k = ngs_sample_sheet.loc[ngs_sample_sheet['plate'] == k]
        # Function to grep barcodes in known positions
        if len(mate1) > 0:
            fastq_sheet[k] = ml.align_fastq_with_plate_positions(mate1, mate2, sample_sheet_plate_k)
# Concatenate dataframes for each plate
fastq_sheet = pd.concat(fastq_sheet)

# Merge file paths with plantIDs and pedigree
sample_sheet = pd.merge(fastq_sheet, pedigree, how = 'left', on = "plantID")
# Explicitly add blank samples as belonging to a blank crossID
sample_sheet.loc[sample_sheet['crossID'].isna(), 'crossID'] = "blank"

# Add paths to the genomes to map to.

# Dictionary lining up possible cross IDs with the paths to the genomes that
# each should be mapped to.
# Note that the direction changes in the keys, but not necessarily in the values.
# See scripts in 03_processing/01_prepare_genomes/ for details of how these were done.
genome_dir = "03_processing/01_prepare_genomes/"
genomes = {
    # Parents
    '1158' : genome_dir + "1158/1158.fa",
    '6024' : genome_dir + "6024/6024.fa",
    '6184' : genome_dir + "6184/6184.fa",
    '8249' : genome_dir + "8249/8249.fa",
      
    '1158x6024' : genome_dir + '1158x6024/1158x6024.fa',
    '6024x1158' : genome_dir + '1158x6024/1158x6024.fa',
    '1158x6184' : genome_dir + '1158x6184/1158x6184.fa',
    '6184x1158' : genome_dir + '1158x6184/1158x6184.fa',
    
    '8249x6024' : genome_dir + '8249x6024/8249x6024.fa',
    '6024x8249' : genome_dir + '8249x6024/8249x6024.fa',
    '8249x6184' : genome_dir + '8249x6184/8249x6184.fa',
    '6184x8249' : genome_dir + '8249x6184/8249x6184.fa',
    # Diagonal crosses
    '1158x8249' : genome_dir + '1158x8249/1158x8249.fa',
    '8249x1158' : genome_dir + '1158x8249/1158x8249.fa',
    '6184x6024' : genome_dir + '6184x6024/6184x6024.fa',
    '6024x6184' : genome_dir + '6184x6024/6184x6024.fa',
    # Blank and controls
    'blank' : '03_processing/01_prepare_genomes/tair10/TAIR10_plus_vectors.fa',
    'Col0'  : '03_processing/01_prepare_genomes/tair10/TAIR10_plus_vectors.fa'
}
# Genomes for each individual
sample_sheet['genome'] = [ genomes[x] for x in sample_sheet['crossID'] ]

# Subset columns and write to disk.
sample_sheet[['plate', 'plantID', 'fastq1', 'fastq2', 'genome']].\
    to_csv('03_processing/04_align_reads/output/sample_sheet.csv', index = False)