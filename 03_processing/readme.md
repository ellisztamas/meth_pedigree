Scripts for processing of raw data, with the output of that code.

* 01_prepare_genomes: Concatentate FASTA files and bisulphite convert
* 02_check_spikeins: Quick checks of sequencing spike-in datasets.
    The output is empty to save disk space
* 03_process_pedigree: Infer ancestors and crossID for each individual
* 04_validate_sample: Align raw bisulphite reads to the TAIR10 genome and
    validate they are who I think they are with ibdpainting.
* 05_genome_assembly: De novo assembly of genomes for accessio s 1158, 6184, 6024
    and 8249 from PacBio HiFi reads.
