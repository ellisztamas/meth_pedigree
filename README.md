# Inheritance of methylation through a pedigree

An experiment looking at large-scale crosses between accessions of *A. thaliana* with high and low TE loads. We anticipate that these data will be used for several projects by separate people, so this repo is intended to document the experimental setup and data processing of reads, which can then be linked to from other project folders.

This README should serve as a quick overview of what we've done. See the [eLabJournal](https://vbc.elabjournal.com/members/experiments/?projectID=14726&studyID=45798) entry for this project has all the notes we made about the experimental procedure.

The github repo is [here](https://github.com/ellisztamas/meth_pedigree), where you can also find an easier-to-read version of this readme file.

## Contents

- [Design](#design)
	- [What to cross](#what-to-cross)
	- [Crosses and sowing](#crosses-and-sowing)
- [Seeds and tissue](#seeds-and-tissue)
- [Data](#data)
	- [Pedigree file](#pedigree-file)
	- [Sequencing data](#sequencing-data)
- [Dependencies](#dependencies)
- [Author information](#author-information)
- [Acknowledgements](#acknowledgements)

## About this folder

Subsfolders are arranged in approximately the order you would need to deal with them:

- `01_raw_data`: soft links to the raw data, which are stored somewhere else
- `02_library`: Directory to store library code, i.e. functions and generalisable scripts that
you would run more than once that you want to apply to multiple datasets, but
keep in one easy-to-maintain place.
- `03_randomisation`: scripts and files related to randomising positions of plants across each experiment
- `04_pedigree`: files related to tracking the family tree of each cross
- `05_bisulphite_data`: Processing of bisulphite reads
- `06_rna_seq`: Processing RNA seq reads

See further readme files within each folder.

## Design

We are interested in what happens when cross accessions with many TEs with accessions with few TEs. As a control, we would also cross high x high and low x low. To try and draw this:

```
High --- High
 |        |
Low ---- Low
```

The idea is to follow methylation status from the parents, to at least the F3 generation. In parallel we also keep track of the selfed offspring of the parents at each generation (henceforth S1, S2, S3; where S stands for selfed).

We have germinated seeds at 21°C for two weeks. Thereafter we grow plants at 10°C, because this better reflects the growing conditions experienced in the field.

### What to cross
 
I (Tom Ellis) picked the following combinations of accessions, plus their reciprocals:

- 6024 x 8249 (low x low TE)
- 6024 x 1158 (low x high TE)
- 8249 x 6184 (low x high TE)
- 1158 x 6184 (high x high TE)

The reason for picking these accessions were that (at the time in 2019 at least):

- They had extreme TE loads. This was based on Mayela's TE calls, which at the time of writing this (Feb 2022) are being revised by Anna Igolkina.
- They are all from the sam region (Sweden)
- Long-read genomes were available
- They have similar flowering times
- Except for DDM1, these crosses do not segregate at any of the major-effect genome-wide loci (NRPE1, MIR823A, various CMT2 alleles identifed by Eriko, CMT3, MET1); we wanted these to not be segregating so that we would have better power to look for smaller-effect loci and changes in epigenetic inheritance.
- There is no evidence for residual heterozygosity

### Crosses and sowing

Crosses were done by Almudena, because Tom and Rahul are clearly much too clumsy to do this. Full details are [here](https://docs.google.com/spreadsheets/d/19ablnFAWTIsy89xEEHPjmm2zuclqxLc2tJ99TpIOBew/edit#gid=0).

We found two issues with cross viability:

1. At the initial cross, crosses involving 1158 did not do well, which Almudena thinks was because she did the emasculation and ferilisation on the same day. She repeated the crosses and it worked much better.
2. Germination for subsequent generations has been patchy. We think that think is because at 10°C seeds take a really long time to dry out, and perhaps we didn't leave some of them long enough.

To address issue 1, Almudena repeated the crosses, so there are two cohorts of crosses.

#### Cohort 1

- Parents:
	- transplanted 5 plants each to pots from seedlings sown by Dejan Dukic for his ATAC-seq project.
	- Two crosses gave very few seeds: 6184x1158 and 1158x6024
	- Seedlings were transplanted 1.7.2019
- F1s:
	- Sowed 6 replicates of each of the 8  (4x2reciprocal crosses) F1 families.
	- Plus 6 replicates of each of 8 S1 families (there were two parents from each parental accession in different crosses)
	- This gives us 6x16=96 pots.
	- Randomisation is done with the R script `F1_randomisation_spring2020.R`.
	- Seeds sown 11.2.2019
- F2s:
	- For each F2 cross we want ~200 individuals.
	- There aren't enough seeds from 6024x1158 or any S2s from 1158, so these aren't included here.
	- Randomistation is in 15 blocks of 96 plants, corresponding to one sequencing plate (48 pots fit in a tray, so two trays is one block)
	- Each block contains 13 individuals from each of F2 7 families (7x13 = 91)
	- The remaining 75 pots (5 per block) are filled with S2s, and their positions are ranomised over the whole experiment.
	- Seeds sown 4.11.2019
	- See `003.scripts/F2_randomisation_november2020.R` for details of the randomisation.
	- For blocks 11 to 15 we discarded plants after collecting tissue, because we did not plan to collect seeds from these (we only wanted them as F2s for mapping, and there are plenty of seeds from the other trays).
- F3s:
	- TBC

#### Cohort 2

- Parents:
	- Sowed plants of all four parental accessions in parallel with the F1s for cohort 1.
	- Almudena repeated the 4 reciprocal crosses, plus the two additional reciprical crosses not done the first time round.
- F1s:
	- Sowed seeds for all 12 reciprocal crosses and 8 S1 families
	- Sown June 2020
	- See [eLabJournal](https://vbc.elabjournal.com/members/experiments/browser/#view=experiment&nodeID=233336&page=0&userID=0&status=0&column=created&order=DESC&search=) for the IDs of families used.

## Seeds and tissue

### Tissue harvesting

- This is somewhat challenging because we need both tissue and seeds from the same plant, cutting off leaves for tissue is likely to affect methylation.
- We harvested tissue as soon as possible after bolting with the idea somatic changes at the rosette are less likely to affect cells in flowers and hence be passed on to offspring. Rahul Pisupati ran a control experiment to check if that is true.
- We took one leaf for bisulphite sequencing and plated directly to sequencing tubes on dry ice. When a row of 8 tubes was full, we dropped this into liquid nitrogen until the all rows in the plate were complete, then moved the plate to -80°C.
- For parents, F1s, and 4 plates of the F2s in cohort 1 we sampled additional leaves for further analyses.

### Additional tissue 

Extra tissue for the parents, F1s and S1s are in boxes in Tom's drawer (row 6, furthest drawer to the right) of -80°C freezer LA00'940'02A.
There is also tissue for four plates of F2s/S2s from blocks 1-4 of cohort 1 in rack C6 (5th row from the bottom, middle row).

### Seeds

F1, and S1 are in the cold in a box labelled 'Methylation pedigree crosses'.
Right now, F2 seeds from cohort 1 are still in bags in the storage room next to the lab waiting to be cleaned.

## Data

### Pedigree file

This project relies on being able to trace the ancestry of each individual, so perfect knowledge of the pedigree is essential.
This is documented in `pedigree.xlsx`, which has been periodically uploaded to [google drive](https://docs.google.com/spreadsheets/d/1jNMV4UUk0v8y47F9XQFN8sAWN-Dn_Rx3IJ5OhCuruHU/edit#gid=0).

The file has two tabs, one listing every individual, and the other summarising families. They show:

- `plantID`: A unique ID for each plant
- `alternative`: For S1s only, and for the family tab, there are two naming systems; see below.
- `mother`: ID of the mother.
- `father`: ID of the father
- `type`: Whether the plant is a parent, from a cross (F1, F2, F3) or selfed offspring of the parents (S1, S2, S3, ...)
- `generation`: Number of generations since the parent. The parents are generation 0.
- `sown`: Date sown, for individuals
- `date_crossed`: Date crossed, for families
- `date_harvested`: Date we collected seed, for families

The `plantID`s need some explanation. Most individuals have a three part name, such as F2.32.169, where:

	- F2 indicates the plant is an F2
	- 32 is the specific family, meaning that all plants labelled F2.32 come from the same mother.
	- 169 is a unique identifier for the individual plant in that family.

Parents are labelled by their accession, and a subscript. These IDs are taken from the notes made by Almudena when she did the crosses. See [here](https://docs.google.com/spreadsheets/d/19ablnFAWTIsy89xEEHPjmm2zuclqxLc2tJ99TpIOBew/edit#gid=0)

For F1s for cohort 1 I initially started with a naming system that I realised wouldn't be feasible as the pedigree got bigger.
Instead they have names like `L2.4 x H2.1 - 5`, where:

- `L2` and `H2` indicate the parents:
	- 6024: L1
	- 8249: L2
	- 1158: H1
	- 6184: H2
- `.4` and `.1` indicate specific individuals of each parent
- `- 5` indicates the specific individual within that F1 family.

### Sequencing data

#### Raw data 

Sequencing plates (see the group NGS [masterlist](https://docs.google.com/spreadsheets/d/1XjO8zabj-1vlu-ex37MRnsnXB_c1U3_k-uXoeKaeGn0/edit#gid=26733257)):

- Parents and F1s are in plate 2021-007
- F2s for cohort 1 are plates 2021-013 to 2021-027

These are being sequenced right now - update this section when we know where the data are.

### Processed sequence data

## Dependencies

None yet - update as necessary.

## Author information

- Experimental design and set-up: Tom Ellis, Rahul Pisupati, Grégoire Bohl Viallefond
- Crosses done by Almudena Morales Mollá
- Principal investigator: Magnus Nordborg

## Acknowledgements

Thanks to Ortun Mittelstein-Schied, Fred Berger and Arturo Marí-Ordonez for useful feedback, and to Joanna Gunis for help with plant work