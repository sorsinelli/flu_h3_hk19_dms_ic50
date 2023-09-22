# Generating targeted mutations for library mutagenesis

Code for generating a file of mutations we want to include in our barcoded DMS libraries, for any given H3 sequence. 

This analysis takes into account: 

* A list of deleterious mutants identified by DMS of Perth09, output from 
[initial_prefs_analysis](https://github.com/jbloomlab/barcoded_H3_DMS/tree/main/library_design/initial_prefs_analysis/results).

* Analysis of whether any of these deleterious mutants are present in any natural H3 sequences circulating since 1968. 
This analysis was conducted by John Huddleston and can be found 
[here](https://github.com/huddlej/barcoded_H3_DMS_natural_frequencies). The alignments included 7,597 variants after basic 
quality control and excluding egg-passaged strains. There were 45 mutations to keep identified, ranging from frequencies 
of 0.5% to 79%. These mutations are listed in `data/mutations_to_keep.csv`. 

* `.fasta` files with nucleotide sequences from Kansas/14/2017 (GISAID accession #WSS2413637) and 
HongKong/45/2019 (GISAID accession #WSS2413591), our two strains of interest.

[h3_chimeric_library_design.ipynb](https://github.com/jbloomlab/barcoded_H3_DMS/blob/main/library_design/targeted_primers/h3_chimeric_library_design.ipynb) 
contains the full analysis.

```python

```
