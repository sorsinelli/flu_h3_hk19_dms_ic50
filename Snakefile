"""Top-level ``snakemake`` file that runs analysis."""


import os


configfile: "config.yaml"


# include `dms-vep-pipeline` pipeline Snakemake file
include: os.path.join(config["pipeline_path"], "pipeline.smk")


rule all:
    input:
        variant_count_files,
        rules.check_adequate_variant_counts.output.passed,
        antibody_escape_files,
        (
            [config["muteffects_observed"], config["muteffects_latent"]]
            if len(func_selections)
            else []
        ),
        config["docs"],


# Arbitrary other rules should be added here
rule get_perth2009_data:
    """Get data from Lee and Eguia studies of Perth/2009 (H3N2) HA."""
    params:
        eguia_url="https://raw.githubusercontent.com/jbloomlab/map_flu_serum_Vietnam_H3_Perth2009/master/results/avgdiffsel/avg_sel_tidy.csv",
        lee_url="https://raw.githubusercontent.com/jbloomlab/map_flu_serum_Perth2009_H3_HA/master/results/avgdiffsel/avg_sel_tidy.csv",
    output:
        eguia_csv="results/perth2009/eguia_avg_sel_tidy.csv",
        lee_csv="results/perth2009/lee_avg_sel_tidy.csv",
    conda:
        "environment.yml"
    log:
        "results/logs/get_perth2009_data.txt",
    shell:
        """
        curl -s {params.eguia_url} > {output.eguia_csv} 2> {log}
        curl -s {params.lee_url} > {output.lee_csv} 2>> {log}
        """


rule process_perth2009_data:
    """Process the Perth/2009 (H3N2) HA data."""
    input:
        rules.get_perth2009_data.output.eguia_csv,
        rules.get_perth2009_data.output.lee_csv,
        nb="notebooks/process_perth2009_data.ipynb",
    output:
        csv="results/perth2009/merged_escape.csv",
        nb="results/notebooks/process_perth2009_data.ipynb",
    conda:
        "environment.yml"
    log:
        "results/logs/process_perth2009_data.txt",
    shell:
        "papermill {input.nb} {output.nb} &> {log}"

   

rule generate_filtered_escape_file:
    """Generate a file of filtered avg escape scores for all HK/19 sera"""
    input:
        expand("results/antibody_escape/{serum}_avg.csv", serum=pd.read_csv(config["antibody_selections"])["antibody"].unique()),
        polyclonal_config=config["polyclonal_config"],
        muteffects_csv="results/muteffects_functional/muteffects_observed.csv",
        nb="notebooks/generate_filtered_escape_file.ipynb",
    output:
        csv="results/full_hk19_escape_scores.csv",
        nb="results/notebooks/generate_filtered_escape_file.ipynb",
    conda:
        "environment.yml"
    log:
        "results/logs/generate_filtered_escape_file.txt",
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


# Add any extra data/results files for docs with name: file
extra_data_files = {
    "sequential to reference site numbering": config["site_numbering_map"],
    "escape for Perth2009 DMS": "results/perth2009/merged_escape.csv",
    "escape for HongKong/2019 DMS": "results/full_hk19_escape_scores.csv",
}


# include `dms-vep-pipeline` docs building Snakemake file
include: os.path.join(config["pipeline_path"], "docs.smk")
