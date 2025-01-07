#!/usr/bin/env bash
# 
# Assemble chloroplast and mitochondria using TIPP.
#
# Note that this uses a conda environment specific for TIPP.
#
# Tom Ellis, 25th November 2024

# SLURM
#SBATCH --job-name=organelle_assembly
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --array=0

# Define working directory and load the conda environment
# (because this is different from what's in setup.sh)
workdir=/scratch-cbe/users/$(whoami)/meth_pedigree
source activate /users/$(whoami)/.conda/envs/TIPP

# === Input === #

i=$SLURM_ARRAY_TASK_ID
# Accession names
accession_array=(1158 6184 8249 6024)
accession_name=${accession_array[$i]}

# Input raw fastq file
input_raw_reads=$workdir/05_genome_assembly/01_raw_data/${accession_name}.fastq.gz

tipp=/groups/nordborg/projects/epiclines/002.pedigree/02_library/tipp_latest.sif

# === Output === #

# Directory to store the output
outdir=$workdir/05_genome_assembly/06_assemble_organelles
mkdir -p $outdir

# Intermediate FASTQ file with downsampled reads
downsampled_reads=$outdir/${accession_name}_downsampled.fastq.gz

# === Script === #

# # Download the container if it doesn't exist already
# cd 02_library
# if [ ! -f tipp_latest.sif ]
# then 
#     mkdir apptainer_cache apptainer_tmp
#     export APPTAINER_CACHEDIR=/groups/nordborg/projects/epiclines/002.pedigree/02_library/apptainer_cache
#     export APPTAINER_TMPDIR=/groups/nordborg/projects/epiclines/002.pedigree/02_library/apptainer_tmp
#     apptainer pull docker://weigelworld/tipp:latest
# fi

cd $outdir

# Downsample the Fastq file so TIPP doesn't melt.
reformat.sh \
    in=$input_raw_reads \
    out=$downsampled_reads \
    samplereadstarget=100k \
    sampleseed=1612


TIPPo.v2.4.pl \
    -t ${SLURM_CPUS_PER_TASK} \
    -g chloroplast \
    -f $downsampled_reads
# singularity exec $tipp /usr/local/bin/TIPP_plastid.v2.1.pl -t ${SLURM_CPUS_PER_TASK} -f $downsampled_reads