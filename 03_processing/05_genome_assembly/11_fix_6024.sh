#!/usr/bin/env bash
# 
# Identify the point of misassembly in the 6024 assembly.
#
# hifiasm creates a large contig (labelled ptg000002l) for accession 6024 that
# spans both chr 2 and 3. I don't think this is real (see 10_linkage_F2_cross.ipynb)
# and I think the contig should be split at position 13865993.

# This script separates contig ptg000002l from all other unscaffolded contigs,
# splits ptg000002l into two, and joins these back to the other contigs.
# It then runs RagTag on the new FASTA file to scaffold the genome.
#
# Tom Ellis, 19th December 2024

# SLURM
#SBATCH --job-name=fix_6024
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --qos=short
#SBATCH --mem=10G
#SBATCH --time=8:00:00

# Define working directory and load the conda environment
source setup.sh


# === Input === #

# FASTA file with unscaffolded contigs
unscaffolded_contigs=$workdir/05_genome_assembly/02_assembly/6024_unscaffolded_contigs.fa

# Complete reference genome to align to
ref_genome=01_data/04_parental_genomes/8334.fa




# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/11_fix_6024
mkdir -p $outdir

# FASTA file with contig ptg000002l only
errant_contig=$outdir/errant_contig.fasta

# Text file with contig names that are *not* ptg000002l
all_other_contig_names=$outdir/ids.txt
# FASTA file with all other contigs
all_other_contigs=$outdir/all_other_contigs.fasta

# SAM file giving alignment positions between the reference genome and the errant contig
alignment_file=$outdir/alignment_file.paf

# FASTA file for separate subcontigs of ptg000002l
subcontig1=$outdir/subcontig1.fasta
subcontig2=$outdir/subcontig2.fasta

# FASTA file containing ptg000002l split into four, plus all other contigs.
new_unscaffolded_contigs=$outdir/new_unscaffolded_contigs.fasta

# FASTA file containing scaffolded contigs matching the autosomes.
scaffolded_contigs=$outdir/scaffolded_contigs.fasta

# === Script === #

# Split the unscaffolded contigs into ptg000002l and everything else.
# Extract the contig from the FASTA of all unscaffolded contigs.
samtools faidx $unscaffolded_contigs ptg000002l > $errant_contig
# Text file listing contigs names *except* ptg000002l
grep ">" $unscaffolded_contigs | grep -v "ptg000002l" | sed 's/^>//' > $all_other_contig_names
# Subset the unscaffolded contigs based on that list.
awk -v ids=$all_other_contig_names 'BEGIN{while((getline<ids)>0)l[">"$1]=1}/^>/{f=l[$1]}f' $unscaffolded_contigs > $all_other_contigs

# Align the errant contig to a complete genome
minimap2 -x asm5 --eqx $ref_genome $errant_contig > $alignment_file

# Split contig ptg000002l
samtools faidx $errant_contig ptg000002l:-13865993 > $subcontig1
samtools faidx $errant_contig ptg000002l:13865994- > $subcontig2
cat $subcontig1 $subcontig2 $all_other_contigs > $new_unscaffolded_contigs


# Scaffold them on down to funky town
ragtag.py scaffold \
    -q 60 \
    -f 10000 \
    -i 0.6 \
    -w \
    --remove-small \
    -o ${outdir} \
    ${ref_genome} \
    $new_unscaffolded_contigs


# Pull out only the scaffolded contigs
samtools faidx \
    $outdir/ragtag.scaffold.fasta \
    Chr1_RagTag Chr2_RagTag Chr3_RagTag Chr4_RagTag Chr5_RagTag > \
    $scaffolded_contigs
samtools faidx $scaffolded_contigs