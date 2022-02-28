Links to the raw data.

Raw reads from the sequencing centre are saved automatically to their own
cluster. In this folder, set up soft links to those folders so we can keep this
project self-contained.

For example, if raw reads for bisulphite reads are save to
`/path/to/bisulphite/data`, set up a soft link with

```
ln -s \
/path/to/bisulphite/data \
01_raw_data/01_bisulphite_reads
```

The image files `cross_data_collection_figure.*` give an overview of the data we planned to collect on different bits of the pedigree as of 11 January 2021. That is likely to have changed by the time you are reading this.