# import basic packages
import pandas as pd
from snakemake.utils import validate


# read sample sheet
local_samples = (
    pd.read_csv(config["local_samples"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

accessions = (
    pd.read_csv(config["accessions"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

samples = local_samples.index.tolist() + accessions.index.tolist()

# validate sample sheet and config file
# validate(samples, schema="../schemas/samples.schema.yaml")
# validate(samples, schema="../schemas/accessions.schema.yaml")
validate(accessions, schema="../schemas/accessions.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

def get_genome_fas(wildcards):
    if wildcards.sample in local_samples.index:
        return local_samples.loc[wildcards.sample, 'assembly_file']
    elif wildcards.sample in accessions.index:
        return f"resources/genomes/{wildcards.sample}.fas"
