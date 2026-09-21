import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")


# report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
configfile: "config/config.yaml"


# read sample sheet
local_samples = (
    pd.read_csv(config["local_samples"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

if local_samples.index.has_duplicates:
    raise ValueError(
        f"Duplicate sample names found in local_samples sheet: {local_samples.index[local_samples.index.duplicated()].tolist()}"
    )

accessions = (
    pd.read_csv(config["accessions"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

if accessions.index.has_duplicates:
    raise ValueError(
        f"Duplicate sample names found in accessions sheet: {accessions.index[accessions.index.duplicated()].tolist()}"
    )


# validate sample sheet and config file
# validate(samples, schema="../schemas/samples.schema.yaml")
# validate(samples, schema="../schemas/accessions.schema.yaml")
validate(local_samples, schema="../schemas/local_samples.schema.yaml")
validate(accessions, schema="../schemas/accessions.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

primers = (
    pd.read_csv(config["primers"], sep="\t", dtype={"Locus": str})
    .set_index(["Locus", "OligoType"], drop=False)
    .sort_index()
)
validate(primers, schema="../schemas/primers.schema.yaml")
if primers.index.has_duplicates:
    raise ValueError(
        f"Duplicate locus name and oligo type combinations found in primers sheet: {primers.index[primers.index.duplicated()].tolist()}"
    )


accessions["assembly_file"] = accessions["sample"].apply(
    lambda sample: f"resources/genomes/{sample}.fas"
)

dataset = pd.concat([local_samples, accessions])
if dataset.index.has_duplicates:
    raise ValueError(
        f"Sample name found in both input sheets: {dataset.index[dataset.index.duplicated()].tolist()}"
    )

# samples = local_samples.index.tolist() + accessions.index.tolist()
samples = dataset.index.tolist()


LOCAL_SAMPLES = "(" + ")|(".join(local_samples.index.tolist()) + ")"
ACCESSION_SAMPLES = "(" + ")|(".join(accessions.index.tolist()) + ")"


def get_genome_fas(wildcards):
    if wildcards.sample in local_samples.index:
        return local_samples.loc[wildcards.sample, "assembly_file"]
    # elif wildcards.sample in accessions.index:
    else:
        return f"resources/genomes/{wildcards.sample}.fas"
