"""Process barcode runs to add samples and ensure runs are unique."""


import pandas as pd


input_csv = snakemake.input.csv
output_csv = snakemake.output.csv
pacbio_libraries = snakemake.params.pacbio_libraries

def process_sample(row):
    if row["library"] not in pacbio_libraries:
        raise ValueError(f"library {row['library']} not in {pacbio_libraries=}")
    label_cols = ["date", "virus_batch", "sample_type"]
    if row["sample_type"] == "antibody":
        label_cols += ["antibody", "antibody_concentration"]
    elif row["sample_type"] not in {"no-antibody_control", "other_control", "plasmid"}:
        raise ValueError(f"invalid `sample_type` {row['sample_type']}")
    label_cols.append("replicate")
    sample = []
    for col in label_cols:
        if pd.isnull(row[col]):
            raise ValueError(f"null {col} in {row}")
        sample.append(row[col])
    return "_".join(map(str, sample))


print(f"Reading barcode runs from {input_csv}")
barcode_runs = (
    pd.read_csv(snakemake.input.csv)
    .assign(
        sample=lambda x: x.apply(process_sample, axis=1),
        library_sample=lambda x: x["library"] + "_" + x["sample"],
    )
)

# ensure sample and FASTQ are unique
for col in ["library_sample", "fastq_R1"]:
    dups = (
        barcode_runs
        .groupby(col)
        .aggregate(n=pd.NamedAgg("library", "count"))
        .query("n > 1")
    )
    if len(dups):
        raise ValueError(f"Found some duplicated {col}:\n{dups}")

print(f"Writing barcode runs with samples to {output_csv}")
barcode_runs.to_csv(output_csv, index=False)
