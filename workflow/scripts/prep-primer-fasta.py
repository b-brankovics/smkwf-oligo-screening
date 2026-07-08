import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w", buffering=1)

input_tsv = snakemake.input["tsv"]
output_fasta = snakemake.output["fasta"]

oligos = pd.read_table(input_tsv, usecols=['Name', 'Sequence'])
with open(output_fasta, "w") as f:
    for index, row in oligos.iterrows():
        f.write(">" + row['Name'] + "\n" + row['Sequence'] + "\n")
