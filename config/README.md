## Workflow overview

This workflow extracts PCR amplicon regions from genomes.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

> Need to update

## Running the workflow

### Input data

This workflow extracts PCR amplicon regions from genomes, and then can process them for STR typing.
You need to specify four tables (TSVs) and a yaml file (primers.yaml) as inputs:

accessions.tsv:

| sample    | assembly        |
| --------- | --------------- |
| Af293	    | GCF_000002655.1 |
| A1160     | GCA_024220425.1 |
| W72310    | GCA_040167795.1 |
| ATCC46645 | GCA_040142955.1 |

local.tsv:

| sample    | assembly_file                      |
| --------- | ---------------------------------- |
| sample1   | data/genomes/sample1_contigs.fasta |
| sample2   | data/genomes/sample2_contigs.fasta |
| sample3   | data/genomes/sample3_contigs.fasta |

primers.tsv:

| Locus    | Description             | Name      | Sequence              | OligoType | Reference                           |
| -------- | ----------------------- | --------- | --------------------- | --------- | ----------------------------------- |
| ITS | Internal Transcribed Spacer | ITS5 | GGAAGTAAAAGTCGTAACAAGG | F | DOI:10.1017/S0953756297005881 |
| ITS | Internal Transcribed Spacer | ITS4 | TCCTCCGCTTATTGATATGC | R | White et al. (1990) |
| LSU | Large subunit ribosomal ribonucleic acid | LROR | ACCCGCTGAACTTAAGC | F | Vilgalys & Hester (1990) |
| LSU | Large subunit ribosomal ribonucleic acid | LR5 | TCCTGAGGGAAACTTCG | R | Vilgalys & Hester (1990) |
| BenA | Beta tubulin | Bt2a | GGTAACCAAATCGGTGCTGCTTTC | F | DOI:10.1128/aem.61.4.1323-1330.1995 |
| BenA | Beta tubulin | Bt2b | ACCCTCAGTGTAGTGACCCTTGGC | R | DOI:10.1128/aem.61.4.1323-1330.1995 |
| CaM | Calmodulin | cmd5 | CCGAGTACAAGGAGGCCTTC | F | DOI:10.1080/15572536.2006.11832738 |
| CaM | Calmodulin | cmd6 | CCGATAGAGGTCATAACGTGG | R | DOI:10.1080/15572536.2006.11832738 |
| Act | nuclear actin | ACT-512F | ATGTGCAAGGCCGGTTTCGC | F | DOI:10.1080/00275514.1999.12061051 |
| Act | nuclear actin | ACT-783R | TACGAGTCCTTCTGGCCCAT | R | DOI:10.1080/00275514.1999.12061051 |
| rodA | hydrophobin | rodA1 | GCTGGCAATGGTGTTGGCAA | F | DOI:10.1080/00275514.1998.12026977 |
| rodA | hydrophobin | rodA2 | AGGGCAATGCAAGGAAGACC | R | DOI:10.1080/00275514.1998.12026977 |
| CYP51A | Erg11 or Cyp51A | cyp51_F | CGGCCGGATGGACATCT | F | DOI:10.1128/jcm.00604-19 |
| CYP51A | Erg11 or Cyp51A | cyp51_R | GCTCGAGCAGCGGTAAAAAT | R | DOI:10.1128/jcm.00604-19 |
| CYP51A | Erg11 or Cyp51A | cyp51_LST | CAATGGCTGAGATTAC | P | DOI:10.1128/jcm.00604-19 |



### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter             | type | details                                             | default                 |
| --------------------- | ---- | --------------------------------------------------- | ----------------------- |
| **accessions**        | path | path to sample sheet of accessions, mandatory       | "config/accessions.tsv" |
| **local_samples**     | path | path to sample sheet of local assemblies, mandatory | "config/local.tsv"      || 