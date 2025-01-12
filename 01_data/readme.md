Overview of raw data files.

These are considered readonly, and any processed data files based on these will
be in 03_processing with the code to generate them.

01_cross_design: Files tracking the initial crosses and subsequent generations
02_grow_plants: scripts to generate randomisation files to grow plants, and the
    files they output.
03_tair10: TAIR10 reference genome
04_parental_genomes: Assemblies for the parental genomes of the crosses
05_spikeins: Data from spikeins for bisulphite runs sorted by date.
06_raw_bisulphite_reads: Raw data for full sequencing runs of bisulphite data,
    sorted by plate. See the Nordborg NGS master list:
    https://docs.google.com/spreadsheets/d/1XjO8zabj-1vlu-ex37MRnsnXB_c1U3_k-uXoeKaeGn0/edit#gid=26733257
07_vector_DNA: Fasta files for the extra vectors used in library preparation
    (phiX, lambda phage, pUC19)
08_1001genomes: Files from the 1001 genomes project, including database files 
    for running SNPmatch using TAIR10 coordinates.
09_pacbio_reads: Raw PacBio HiFi reads for the parental accessions.
10_Col-CEN: Data for telomere to telomere assembly of Columbia from
    Naish et al. Science374,eabi7489(2021).DOI:10.1126/science.abi7489
    To clone:
    git clone https://github.com/schatzlab/Col-CEN.git
