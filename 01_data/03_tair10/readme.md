Files relating to the the TAIR10 genome assembly.

Most of these files are from the official TAIR website:
https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release

`TAIR10.hard_masked.fa` is the version of the TAIR10 genome used by Rabanal et al
for HiFi assemblies which is masked for organellar and repetitive sequences.
See https://doi.org/10.5281/zenodo.7326462

`ColCC.fa` is a pre-release draft of the TAIR12 genome.
This file has unusual FASTA headers; this code changes it to Chr1, Chr2 etc,
and creates `TAIR12.fa`
```
cat 01_data/03_tair10/ColCC.fa | awk '{if(substr($1,1,1) == ">") {printf ">Chr%s\n", $NF;} else {print}}' > 01_data/03_tair10/TAIR12.fa
```
