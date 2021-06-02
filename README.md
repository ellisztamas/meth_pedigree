# Inheritance of methylation through a pedigree

We crossed accessions of *Arabidopsis thaliana* that differ in transposon load and track changes in cytosine methylation through to the F~3~ generation. We can use these data to address several questions:

1. How hertitable (or not) are methylation marks?
2. How does this correlate with chromatin marks and gene expression?
3. Using the (large sample of ) F~2~ plants, what trans-acting modifiers of genome-wide or local methylation are segregating?
4. How frequent are paramutations? i.e. how often do loci that were methylated in one chromosome but not at the same position in the homologous chromosome in the F~2~ become methylated in the F~3~.
5. What happens with transposons when you cross accessions with wildly differing TE loads?

See also the github repository and [eLabJournal](https://vbc.elabjournal.com/members/experiments/?projectID=14726&studyID=45798) entry for this project.

## Contents

1. [Design](#design)
	1. [Experimental design](#experimental-design)
	2. [Analysis workflow](#analysis-workflow)
3.  [Data](#data)
	1. [Raw data](#raw-data)
	2. [Derived data](#derived-data)
4. [Dependencies](#Dependencies)
5. [License](#license)
6. [Author](#author)
7. [Support](#support)
8. [Acknowledgements](#acknowledgements)

## Design
### Experimental design

* Give enough detail of what was done so that somebody could reproduce this.
* Describe what phenotypes were measured. If it is more appropriate, point to a README file for individual datasets.
* Point to files describing randomisation.

### Sample processing

* What did you do with samples after harvesting?
* Where are they stored?

## Data
### Raw data

* Raw data refers to data collected from experiments or raw sequencing reads.
* These data should be stored in the `/groups/nordborg/raw.data` and should be read only, even to you.
* You can make the system behave as if your data was in your project foler by creating a 'soft link' from your project folder to the raw data file. For example:
```
ln -s /groups/nordborg/raw.data/athaliana/dna/my_raw_reads.bam /001.data/001.raw/my_raw_reads.bam
```

### Derived data

* Describe what pipelines or scripts were run on the raw data, where the scripts are to run them, and where the output is saved.

## Dependencies

List the R/python/other packages needed to run the analyses, along with version numbers. If you have used an environment file, mention that.

## License

For example:
> Released under the MIT license. See LICENSE for full details.

List any information about embargoes on the data.

## Author information

This repository: Tom Ellis (thomas.ellis@gmi.oeaw.ac.at)

Generally, give the names of those involved in data collection, processing, analysis and supervision.

## Acknowledgements

Thank anyone who is not an author, but deserves mentioning.