import sys
import pandas as pd
import yaml

sys.stderr = open(snakemake.log[0], "w", buffering=1)

input_tsv = snakemake.input["tsv"]
output_yaml = snakemake.output["yaml"]

oligos = pd.read_table(input_tsv, usecols=['Locus', 'Name', 'OligoType'])
oligo_dict = {
    locus: {
        row['OligoType']: row['Name']
        for index, row in oligos[oligos['Locus'] == locus].iterrows()
    }
    for locus in oligos['Locus'].unique()
}
with open(output_yaml, 'w') as outfile:
    yaml.dump(oligo_dict, outfile, default_flow_style=False)