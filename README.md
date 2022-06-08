# Deep mutational scanning of H3N2 flu strain A/HongKong/45/2019

For detailed and easy-to-read documentation of the workflow and analysis, see [https://jbloomlab.github.io/flu_h3_hk19_dms](https://jbloomlab.github.io/flu_h3_hk19_dms/).

## Running the analysis
The analysis uses the `conda` environment in [environment.yml](environment.yml).
Build this environment and activate it before running the analysis.

The analysis is run by the `snakemake` file in [Snakefile](Snakefile).
To run the analysis using `<NCPUS>` cores, activate the `conda` environment and run:

    snakemake -j <NCPUS>

To run on the Hutch cluster on 32 cores, use the script [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch -c 32 run_Hutch_cluster.bash 

## Configuring the analysis
The configuration for the analysis is in [config.yaml](config.yaml).
The code that is run by [Snakefile](Snakefile) consists of Jupyter notebooks in [./notebooks/](notebooks) and Python scripts in [./scripts/](scripts).
Input data are in [./data/](data).
Results are placed in [./results/](results) are are **not** tracked by `git` except for the specific exceptions stipulated in [.gitignore](.gitignore).

The pipeline in [Snakefile](Snakefile) also builds HTML documentation of the analysis, which is placed in [./docs/](docs) and can be viewed via GitHub Pages.
See [Rendering docs](#rendering-docs) for more detail.

The Jupyter notebooks in [./notebooks/](notebooks) should have their output stripped by [nbstripout](https://github.com/kynan/nbstripout) before being committed with `git`.
The easiest way to do this is installing [nbstripout](https://github.com/kynan/nbstripout) to your current copy of the repo with:

    nbstripout --install

## Rendering docs
The documentation for the analysis is built by rendering Jupyter notebooks to HTML and including links to key data files.
The final few rules in [Snakefile](Snakefile) build a [sphinx](https://www.sphinx-doc.org/) index file ([./docs_source/index.rst](docs_source/index.rst)) that is then rendered to HTML documentation in [./docs/](docs), which is displayed [via GitHub pages](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site) at [https://jbloomlab.github.io/SARS-CoV-2_Delta_spike_DMS](https://jbloomlab.github.io/SARS-CoV-2_Delta_spike_DMS).
This requires [sphinx](https://www.sphinx-doc.org/), [nbsphinx](https://nbsphinx.readthedocs.io/), and [nbpshinx-link](https://nbsphinx-link.readthedocs.io), all of which are installed in the `conda` environment in [environment.yml](environment.yml).

The relevant rules at the end of [Snakefile](Snakefile) build HTML versions of all Jupyter notebook outputs of prior rules that are named `nb`, and link them in [./docs_source/index.rst](docs_source/index.rst).
They also build and link a rule graph for the `snakemake` workflow, and include links to the data files specified in the `data_files` variable in [Snakefile](Snakefile).
Note these data files must also be specified for inclusion in the GitHub repo via [.gitignore](gitignore).

The directory [./docs_source](docs_source) also contains the other [sphinx](https://www.sphinx-doc.org/) configuration files for building the documentation (i.e., [./docs_source/conf.py](docs_source/conf.py)), and [Makefile](Makefile) is the [sphinx]((https://www.sphinx-doc.org/) file that builds the documentation in [./docs_source/](docs_source) into the HTML.
These files are manually edited versions of ones originally created by the [sphinx]((https://www.sphinx-doc.org/)) `quickstart` command.

## `git lfs` for large files
A few of the large results files are tracked using [git lfs](http://arfc.github.io/manual/guides/git-lfs): e.g., [results/barcodes_runs/barcode_counts.csv](results/barcode_runs/barcode_counts.csv).
Be sure to also track any other large files by `git lfs` rather than regular `git`.

## Library design
The code and data used to design the mutant library are in [./library_design/](library_design), and are run outside the main pipeline in [Snakefile](Snakefile).
