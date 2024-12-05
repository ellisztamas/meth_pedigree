#!/usr/bin/env bash
# 
# Scaffold raw assemblies to:
#     TAIR10 (hard masked to remove repetive sequences)
#     TAIR12 (draft)
#     The telomere-to-telomere assembly from Naish et al (2021)
#     CLR assemblies for the same genome
# After scaffolding, only scaffolds matching autosomes are retained.
# This does not work for the CLR assemblies, because the chromosome names are
# different.
#
# Input:
#     Unscaffolded FASTA file from hifiasm, derived from the 'main' GFA file
#         (ending in .bp.p_ctg.gfa) created by hifiasm
#     FASTA files for the four genomes to scaffold to.
# Output:
#     Output files from ragtag
#     FASTA file of scaffolded contigs for the autosomes only.
#       
# Tom Ellis, 2nd December 2024, adapting code by Robin Burns.

# SLURM
#SBATCH --job-name=scaffolding
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=rapid
#SBATCH --time=30:00
#SBATCH --array=0-3

date

# Define working directory and load the conda environment
source setup.sh


# === Input === #

i=$SLURM_ARRAY_TASK_ID

# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Input raw bam file
unscaffolded_contigs=$workdir/05_genome_assembly/02_assembly/${accession_name}_unscaffolded_contigs.fa
echo "Scaffolding contigs from: ${unscaffolded_contigs}"

# Reference genomes to scaffold to
tair10=01_data/03_tair10/TAIR10.hard_masked.fa
tair12=01_data/03_tair10/TAIR12.fa
col_cen=01_data/10_col-CEN/v1.2/Col-CEN_v1.2.fasta.gz
lu1=01_data/04_parental_genomes/1002.fa
ale_stenar=01_data/04_parental_genomes/8334.fa
# Find the CLR assembly. Luckily, the files are in the same order as the array above.
clr_assembly_array=(01_data/04_parental_genomes/*fasta)
clr_assembly=${clr_assembly_array[$i]}
# Put the paths in an array we can loop over
ref_genome_array=($tair10 $tair12 $col_cen $clr_assembly $lu1 $ale_stenar)


# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/03_scaffolding/$accession_name
mkdir -p $outdir

# RagTag will only take a dircetory name as an output prefix, so define 
# separate output subdirectories for each.
# This is an array, because there's a loop later
output_dir_array=('tair10' 'tair12' 'col_cen' 'clr_assembly' '1002' '8334')

# There's another output file for scaffolded contigs, but the path depends on 
# the arrays above, so it was easier to define in inside the loop.



# === Script === #

# Scaffolding raw contigs for a single accession is done as a loop over four reference genomes.
# It creates separate directories for each reference genome.
for a in {1..5}
do 
    # Define input and outputs
    ref_genome=${ref_genome_array[$a]}
    output_prefix=${outdir}/${output_dir_array[$a]}
    scaffolded_contigs=$output_prefix/${accession_name}_${output_dir_array[$a]}_scaffolded_contigs.fasta
    echo "Scaffolding to reference genome ${ref_genome}"
    echo "Saving to ${output_prefix}"

    # If you have run RagTag before and there are files left over, RagTag seems to
    # use these, even if you use the -w flag to overwrite.
    # To be on the safe side, remove the whole directory.
    if [ -d $output_prefix ]; then
        echo "Removing existing output directory to avoid RagTag caching things."
        rm -r $output_prefix
    fi
    
    # Scaffolding
    SECONDS=0
    ragtag.py scaffold \
        -q 60 \
        -f 10000 \
        -i 0.6 \
        -w \
        --remove-small \
        -o ${output_prefix} \
        ${ref_genome} \
        $unscaffolded_contigs
    duration=$SECONDS
    echo -e "Finished in $((duration / 60)) minutes and $((duration % 60)) seconds."

    # Pull out only the scaffolded contigs
    echo "Pulling out only the contigs matching autosomes, and saving as: ${scaffolded_contigs}"
    samtools faidx \
        $output_prefix/ragtag.scaffold.fasta \
        Chr1_RagTag Chr2_RagTag Chr3_RagTag Chr4_RagTag Chr5_RagTag > \
        $scaffolded_contigs
    samtools faidx $scaffolded_contigs
    echo -e "Finished."
    echo ""
    echo ""

done

date