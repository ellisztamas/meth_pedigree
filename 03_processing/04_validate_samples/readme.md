Script to check that the samples are what I think they are.

The general idea is to align raw reads to the TAIR10 reference genome, call
genotypes at SNPs that vary between the four parents, then to plot how similar
each sample is to what I think the parents are.

An incomplate list of issues:
- 6184_2 is given twice in plate 2021-007 positions D7 and E7.
    The second of these has poor DNA.
    I suspect that the second was meant to be 6184_4, and it seems I realised 
    this was missing and added it later in well D10.
    I will use the data from D7 for 6184_2 for now.
    Consider including this in a later run!
- S1.07.03 is given twice in plate 2021-007 positions H10 and E12.
    I don't know why.
    It doesn't have any offspring, so no problem.

- Well given as empty in the NGS master list:
    - F2.27.050
    - F2.27.052
- Fastq file empty or missing:
    - F3.18.003
    - F3.21.002