# Generating primers for targeted library mutations

This dir contains the collected inputs for running `create_primers.py` in [TargetedTilingPrimers](https://github.com/jbloomlab/TargetedTilingPrimers/), as well as analysis of the resulting set of primers.

[data](https://github.com/jbloomlab/barcoded_H3_DMS/tree/main/library_design/targeted_library_primers/data):

* Codon frequency table for Influenza H3N2, pulled from [kazusa.or.jp Codon Usage Database](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=41857). This is used for codon optimization.

* `fasta` files of the *chimeric* coding sequences for Kansas/2017 and HongKong/2019. These sequences are plasmids [3022](https://github.com/jbloomlab/plasmids/blob/master/genbank_maps/3022_pHH_WSNHAflank_KS17-cterm-nopac_WSNHAduppac-stop_edited.gb) and [3023](https://github.com/jbloomlab/plasmids/blob/master/genbank_maps/3022_pHH_WSNHAflank_KS17-cterm-nopac_WSNHAduppac-stop_edited.gb), respectively. Full details of the chimeric structures can be found in the descriptions of these genbank files. The wildtype ectodomain sequence is maintained apart from the first two amino acids, where a poly-A run has been edited out. 

* `csv` files containing a set of targeted library mutations for Kansas/2017 and HongKong/2019.

`primer_input_data_prep.ipynb`: formats this input data according to the specifications of the TargetedTilingPrimers script, and prints the formatted files to [results/TargetedTilingPrimers_inputs](results/TargetedTilingPrimers_inputs). 

[results](https://github.com/jbloomlab/barcoded_H3_DMS/tree/main/library_design/targeted_library_primers/results):

* `csv` files containing the full set of primer sequences for Kansas/2017 and HongKong/2019.