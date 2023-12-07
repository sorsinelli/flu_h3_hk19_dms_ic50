# Deep mutational scanning of H3N2 flu strain A/HongKong/45/2019

For detailed and easy-to-read documentation of the workflow and analysis, see [https://dms-vep.github.io/flu_h3_hk19_dms](https://dms-vep.github.io/flu_h3_hk19_dms/).

## Organization of this repo

### Summary escape files

A csv file containing the fully processed, filtered escape data is available at [results/full_hk19_escape_scores.csv](results/full_hk19_escape_scores.csv). If you'd like to use the DMS escape data for any additional analyses, this file is the most convenient format. Data has been filtered to only include mutations that were seen in at least 3 distinct variants, are not highly deleterious, and are not stop codons. Columns are as follows - 
* *site*: site numbering in quasi-H3 HA reference numbering. We use the standard H3 HA numbering scheme, where the first codon of the ectodomain is site 1, through the HA1 subdomain. Rather than switching to (HA2)1 at the beginning of the HA2 domain, we continue with sequential numbering (so site (HA2)1 in reference numbering is site 330 in this data). 
* *wildtype*, *mutant*, *mutation*: wildtype amino acid and relevant mutation at the given site.
* *escape_mean*: escape score calculated for that mutation, averaged between two independent libraries. Positive values indicate that the mutation confers escape from the given serum. Negative values indicate that the mutation makes the virus more sensitive to neutralization by that serum.
* *escape_std*: the standard deviation in escape scores between the two libraries.
* *n_models*, *frac_models*: number/fraction of libraries that mutation was represented in. The two libraries were independently mutagenized, so some mutations are only seen in a single library, although most are seen in both libraries.
* *times_seen*: number of unique variants (summed across both libraries) that express the given mutation. This file only includes mutations with `times_seen>=3`
* *functional_effect*: effect of the mutation on HA protein function, calculated based on the ratio of mutation counts in the plasmid library vs the passaged virus library. This file only includes mutations with `functional_effect>=-1.38`. 
* *serum*: numerical ID of the serum sample for the given escape scores.
* *cohort*: age cohort for that serum sample - either `2-5_years` (8 sera), `15-20_years` (8 sera), or `40-45_years` (10 sera). 

This repository also includes finalized analyses of serum escape from a non-barcoded library in the background of A/Perth/16/2009. Because of the difference in library design, the bulk of the data analysis can be found in the separate repository [map_flu_serum_Vietnam_H3_Perth2009](https://github.com/jbloomlab/map_flu_serum_Vietnam_H3_Perth2009/tree/master). However, the finalized, fully processed escape data is available here at [results/perth2009/merged_escape.csv](results/perth2009/merged_escape.csv). Columns are as follows - 
* *name*: serum sample for the given escape scores, given as age and country where serum was collected. This file includes some miscellaneous adult sera collected in Seattle, but all downstream analysis only uses the Vietnam adult and child sera, plus the ferret sera.
* *serum*: numeric ID for that serum.
* *serum_group*: equivalent to *cohort* for HK/19 data - either `adult`, `child`, `VIDD sera` (for misc. Seattle adults), or `ferret`.
* *epitope*: escape is calculated using the software package [polyclonal](https://jbloomlab.github.io/polyclonal/index.html), which attempts to fit *n* epitopes to the data. Here, we specify that `n_epitopes=1`, so this column can be ignored.
* *site*: site numbering in H3 HA reference numbering, where sites in the HA2 subdomain are numbered as `(HA2)1, (HA2)2...` etc.
* *site_sequential*: sites numbered sequentially such that the start codon is 1, as opposed to -16 in H3 HA reference numbering.
* *wildtype*, *mutant*: wildtype amino acid and relevant mutation at the given site.
* *escape*: median escape score calculated for that mutation from three independent libraries. Note that the magnitude of escape scores for the Perth/2009 are not directly comparable between sera, as the experiments are not normalized against a neutralization standard. These experiments are also unable to resolve neutralization-sensitizing mutations, i.e. negative escape scores. See [map_flu_serum_Vietnam_H3_Perth2009](https://github.com/jbloomlab/map_flu_serum_Vietnam_H3_Perth2009/tree/master) for detailed analysis generating these escape scores.

### `dms-vep-pipeline` submodule

Most of the analysis is done by the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline](dms-vep-pipeline), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline` submodule.
You can also make changes to the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) that you commit back to that repo.

### Code and configuration
The [snakemake](https://snakemake.readthedocs.io/) pipeline itself is run by the [Snakefile](Snakefile), which includes [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) reads its configuration from [config.yaml](config.yaml).
The [conda](https://docs.conda.io/) environment used by the pipeline is that specified in the `environment.yml` file in [dms-vep-pipeline](dms-vep-pipeline).

Additional scripts and notebooks that are specific to this analysis and not part of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) are in [scripts/](scripts) and [notebooks/](notebooks).

### Input data
Input data for the pipeline are in [data/](data).

### Results and documentation
The results of running the pipeline are placed in [results/](results).
Only some of these results are tracked to save space (see [.gitignore](.gitignore)).

The pipeline builds HTML documentation for the pipeline in [docs/](docs), which is rendered via GitHub Pages at [https://dms-vep.github.io/flu_h3_hk19_dms](https://dms-vep.github.io/flu_h3_hk19_dms/).

### Library design
The design of the mutant library is contained in [library_design/](library_design).
That design is not part of the pipeline but contains code that must be run separately with its own [conda](https://docs.conda.io/) environment.

## Running the pipeline
To run the pipeline, build the conda environment `dms-vep-pipeline` in the `environment.yml` file of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), activate it, and run [snakemake](https://snakemake.readthedocs.io/), such as:

    conda activate dms-vep-pipeline
    snakemake -j 32 --use-conda --rerun-incomplete

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch -c 32 run_Hutch_cluster.bash
