"""Implements ``snakemake`` rule `run_nb.py` to run Jupyter notebook in place."""


import os
import subprocess


input_nb = snakemake.input.nb
output_nb = snakemake.output.nb

print(f"Running notebook {input_nb} to create {output_nb}")

if os.path.basename(input_nb) != os.path.basename(output_nb):
    raise ValueError(f"{input_nb=} and {output_nb} have different base names")

subprocess.check_call(['jupyter', 'nbconvert',
                       '--to', 'notebook',
                       '--execute',
                       '--output', os.path.abspath(output_nb),
                       '--ExecutePreprocessor.timeout=-1',
                       input_nb,
                       ])