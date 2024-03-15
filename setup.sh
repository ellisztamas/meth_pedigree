#!/usr/bin/env bash

# Assign variables common to all data processing scripts.
# If you need to change something, change it here.

# Tom Ellis, 31st January

echo "Loading working directory and Conda environment from setup.sh"

# Data processing was done on a special drive for large jobs on the VBC CLIP cluster
# This won't exist on other machines, so change it here if you have downloaded
# the code somewhere else.
workdir=/scratch-cbe/users/$(whoami)/meth_pedigree

# workdir=03_processing/pieters_sample_sheet/output # example alternative working directory
echo "Working directory: ${workdir}" 
mkdir $workdir -p

# Load conda environment
# This is probably also CLIP-specific
# If you haven't already, install the environment with `conda env create -f environment.yml`
source activate /users/$(whoami)/.conda/envs/epiclines