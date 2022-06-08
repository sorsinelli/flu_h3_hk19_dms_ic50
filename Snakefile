"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os
import textwrap

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: "config.yaml"

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config["pacbio_runs"], dtype = str)
               .assign(pacbioRun=lambda x: x["library"] + "_" + x["run"].astype(str))
               .set_index("pacbioRun")
               )
assert pacbio_runs.index.nunique() == len(pacbio_runs)

# Rules ---------------------------------------------------------------------
rule all:
    input:
        "docs",

rule gene_sequence:
    """Get the sequence of the gene used in the experiments."""
    input: **config["amplicons"]
    output:
        codon=config["gene_sequence_codon"],
        prot=config["gene_sequence_protein"]
    conda: "environment.yml"
    script: "scripts/gene_sequence.py"

rule align_parse_PacBio_ccs:
    """Align and parse PacBio CCS FASTQ file."""
    input:
        fastq=lambda wc: pacbio_runs.at[wc.pacbioRun, "fastq"],
        amplicon=lambda wc: config["amplicons"][pacbio_runs.at[wc.pacbioRun, "amplicon"]],
        specs=lambda wc: config["amplicon_specs"][pacbio_runs.at[wc.pacbioRun, "amplicon"]],
    output: outdir=directory(os.path.join(config["process_ccs_dir"], "{pacbioRun}"))
    conda: "environment.yml"
    script: "scripts/align_parse_PacBio_ccs.py"

rule analyze_pacbio_ccs:
    """Analyze PacBio CCSs and get ones that align to amplicons of interest."""
    input:
        expand(rules.align_parse_PacBio_ccs.output.outdir,
               pacbioRun=pacbio_runs.index),
        config["amplicons"].values(),
        config["amplicon_specs"].values(),
        nb="notebooks/analyze_pacbio_ccs.ipynb"
    output:
        config["aligned_ccs_files"].values(),
        nb="results/notebooks/analyze_pacbio_ccs.ipynb"
    conda: "environment.yml"
    script: "scripts/run_nb.py"

rule build_pacbio_consensus:
    """Build PacBio consensus sequences for barcodes."""
    input:
        config["aligned_ccs_files"].values(),
        config["gene_sequence_codon"],
        nb="notebooks/build_pacbio_consensus.ipynb"
    output:
        config["nt_variants"].values(),
        nb="results/notebooks/build_pacbio_consensus.ipynb"
    conda: "environment.yml"
    script: "scripts/run_nb.py"

rule analyze_codon_variants:
    """Analyze mutations in variants."""
    input:
        config["nt_variants"].values(),
        config["gene_sequence_codon"],
        config["designed_mutations"],
        config["site_numbering_map"],
        nb="notebooks/analyze_codon_variants.ipynb"
    output:
        config["codon_variants"].values(),
        nb="results/notebooks/analyze_codon_variants.ipynb"
    conda: "environment.yml"
    script: "scripts/run_nb.py"

rule virus_variants_w_neut_standard:
    input:
        config["neut_standard_barcodes"],
        config["codon_variants"]["plasmid"],
        nb="notebooks/virus_variants_w_neut_standard.ipynb"
    output:
        config["virus_variants_w_neut_standard"],
        nb="results/notebooks/virus_variants_w_neut_standard.ipynb"
    conda: "environment.yml"
    script: "scripts/run_nb.py"

checkpoint process_barcode_runs:
   """Check barcode runs and write CSV with sample names and other information."""
   input: csv=config["barcode_runs_wo_sample"]
   output: csv=config["barcode_runs"]
   params: pacbio_libraries=pacbio_runs["library"].tolist()
   conda: "environment.yml"
   script: "scripts/process_barcode_runs.py"

def barcode_count_fate_csvs(wildcards):
   """Get list of all barcode count CSVs."""
   fname = checkpoints.process_barcode_runs.get(**wildcards).output.csv
   library_samples = pd.read_csv(fname)["library_sample"]
   return [
       os.path.join(config[f"barcode_{ftype}_dir"], f + ".csv")
       for f in library_samples for ftype in ["counts", "counts_invalid", "fates"]
   ]

def barcode_fastq_R1(wildcards):
   """Get R1 FASTQ for a specific library-sample."""
   fname = checkpoints.process_barcode_runs.get(**wildcards).output.csv
   return (
       pd.read_csv(fname)
       .set_index("library_sample")
       .at[wildcards.library_sample, "fastq_R1"]
       .split(";")
   )

rule count_barcodes:
   """Count barcodes for a specific library-sample."""
   input:
       fastq_R1=barcode_fastq_R1,
       variants=config["virus_variants_w_neut_standard"],
   output:
       counts=os.path.join(config["barcode_counts_dir"], "{library_sample}.csv"),
       counts_invalid=os.path.join(config["barcode_counts_invalid_dir"], "{library_sample}.csv"),
       fates=os.path.join(config["barcode_fates_dir"], "{library_sample}.csv"),
   conda: "environment.yml"
   script: "scripts/count_barcodes.py"

rule analyze_barcode_counts:
   """Aggregate barcode counts for all samples and make summary plots."""
   input:
       barcode_count_fate_csvs,
       config["barcode_runs"],
       config["gene_sequence_codon"],
       config["neut_standard_barcodes"],
       config["virus_variants_w_neut_standard"],
       nb="notebooks/analyze_barcode_counts.ipynb"
   output:
       config["barcode_counts"],
       nb="results/notebooks/analyze_barcode_counts.ipynb"
   conda: "environment.yml"
   script: "scripts/run_nb.py"

rule analyze_variant_counts:
   """Analyze counts at the variant level for all samples."""
   input:
       config["barcode_counts"],
       config["virus_variants_w_neut_standard"],
       config["gene_sequence_codon"],
       config["barcode_runs"],
       site_map=config["site_numbering_map"],
       nb="notebooks/analyze_variant_counts.ipynb"
   output:
       config["codon_variant_table_pickle"],
       nb="results/notebooks/analyze_variant_counts.ipynb"
   conda: "environment.yml"
   script: "scripts/run_nb.py"

rule variant_prob_escapes:
    """Compute probability of escape for variants."""
    input:
        config["codon_variant_table_pickle"],
        config["site_numbering_map"],
        config["barcode_runs"],
        nb="notebooks/variant_prob_escapes.ipynb"
    output:
        config["variant_prob_escapes"],
        config["antibody_selections"],
        nb="results/notebooks/variant_prob_escapes.ipynb"
    conda: "environment.yml"
    script: "scripts/run_nb.py"

# -------------------------------------------------------------------------
# The following rules must be last in order to properly build nbsphinx docs
# of the `nb` output of all rules that have that output.
# -------------------------------------------------------------------------

# all rules with `nb` as an output are shown as analysis notebooks
nb_output_rules = [name for name, obj in rules.__dict__.items()
                   if hasattr(obj.output, "nb")]

# data files (should also be in `.gitignore`)
data_files = {
    "parental gene sequence (codon)": rules.gene_sequence.output.codon,
    "parental gene sequence (protein)": rules.gene_sequence.output.prot,
    "sequential to reference site numbers": config["site_numbering_map"],
    **{f"{name} barcode to codon-variant table": table
       for name, table in config["codon_variants"].items()},
    "barcode counts": config["barcode_counts"],
    "sample pairings for antibody selections": config["antibody_selections"],
    "variant prob escapes for antibody selections": config["variant_prob_escapes"],
    }

rule make_nblinks:
    """Make `nblink` files for building ``sphinx`` HTML docs of analysis."""
    input: [getattr(getattr(rules, name).output, "nb") for name in nb_output_rules]
    output: [f"docs_source/{name}.nblink" for name in nb_output_rules]
    run:
        for nb, nblink in zip(input, output):
            with open(nblink, "w") as f:
                f.write(f'{{\n    "path": "../{nb}"\n}}')

rule make_rulegraph:
    """Build the ``snakemake`` rule graph."""
    input: "Snakefile"
    output: svg="results/rulegraph.svg"
    conda: "environment.yml"
    shell: "snakemake --rulegraph | dot -Tsvg > {output.svg}"

rule sphinx_index:
    """Create ``index.rst`` file for ``sphinx`` to render HTML docs of analysis."""
    input:
        data_files.values(),
        nblinks=rules.make_nblinks.output,
        rulegraph=rules.make_rulegraph.output.svg,
    output: index="docs_source/index.rst",
    run:
        github_url=f"https://github.com/{config['github_user']}/{config['github_repo']}"
        blob_path=f"{github_url}/blob/{config['github_branch']}"
        indent_space = "\n" + " " * 19  # separate indented entries in tables of contents
        analysis_nbs = indent_space.join([os.path.splitext(os.path.basename(nblink))[0]
                                          for nblink in input.nblinks])
        no_indent_space = "\n" + " " * 16  # separate non-idented entries
        data_file_links = no_indent_space.join([f"- `{label} <{blob_path}/{link}>`_"
                                                for label, link in data_files.items()])
        with open(output.index, "w") as f_obj:
            f_obj.write(textwrap.dedent(f"""\
                {config["description"]}
                {"=" * len(config["description"])}
                Documentation of data analysis.

                For actual code, see {github_url}.
                
                Study by {config["authors"]}.

                Workflow
                --------
                Rule graph for the `snakemake <https://snakemake.readthedocs.io/>`_
                workflow in `Snakefile <{blob_path}/Snakefile>`_:

                .. image:: ../{input.rulegraph}

                Analysis notebooks
                ------------------
                Many of the plots in these notebooks are interactive, so try mousing
                over points for details or zooming.

                .. toctree::
                   :maxdepth: 1

                   {analysis_nbs}

                Key data files
                --------------

                {data_file_links}
                """))

rule sphinx_docs:
    """Build ``sphinx`` docs for project."""
    input:
        rules.sphinx_index.output.index,
        rules.sphinx_index.input,
    output:
        temp_docs=temp(directory("results/docs/html")),
        final_docs=directory("docs"),
    conda: "environment.yml"
    shell:
        """
        make html
        cp -r {output.temp_docs} {output.final_docs}
        touch docs/.nojekyll
        """
