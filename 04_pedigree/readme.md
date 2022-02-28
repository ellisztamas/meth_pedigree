# Pedigree

This folder contains files to keep track of the pedigree (the family tree between all the samples in the project).

## Crosses files

Two files, `191011_parental_crosses.xlsx` and `200610_parental_crosses.xlsx`, are hard copies of the notes about crosses done by Almudena, which she kept on [google drive](https://docs.google.com/spreadsheets/d/19ablnFAWTIsy89xEEHPjmm2zuclqxLc2tJ99TpIOBew/edit#gid=0).

## Pedigree file

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