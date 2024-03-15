Symlinks to fasta files for the genomes of the four parents of the crosses:

- 6024
- 8249
- 6184
- 1158

```
# Code to set up the symlinks
outdir=01_data/04_parental_genomes/
ln -s /groups/nordborg/common_data/eracaps/assemblies/v02.1/AT_6024.v2.1.addChrMC.fasta $outdir
ln -s /groups/nordborg/common_data/pacbio_genomes/assemblies/v01.1/8249_v1.1.fasta $outdir
ln -s /groups/nordborg/common_data/pacbio_genomes/assemblies/v01.1/1158_v1.1.fasta $outdir
ln -s /groups/nordborg/common_data/pacbio_genomes/assemblies/v01.1/6184_v1.1.fasta $outdir
```