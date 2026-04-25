# smkwf-oligo-screening

Snakemake workflow for screening olignucleotide anealing sites and identifying PCR amplicon regions

## TODO

- [ ] Read samples from config files
- [ ] Add schema for checking configs
- [ ] Option to combine local and remote (NCBI) data
- [ ] Use [https://github.com/snakemake-workflows/snakemake-workflow-template](https://github.com/snakemake-workflows/snakemake-workflow-template)
- [ ] Add max length for amplicons
- [ ] Add report as overview for probes and amplicons accross the samples 
- [x] Update for Perl script integration

Bug in input function when using as a module. Could it be that the function doesn't get the same "external" variables? Should it be passed explicitly, or should it use a function to get the same value?
