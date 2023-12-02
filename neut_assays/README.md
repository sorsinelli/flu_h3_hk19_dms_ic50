# Neutralization Assays

This directory includes neutralization data and analysis for both A/HongKong/45/2019 and A/Perth/16/2009. The neutralization data is pre-formatted as fraction infectivity at each serum concentration, which is calculated by normalizing GFP fluorescence in serum+virus+cells wells against GFP fluorescence in control virus+cells wells (i.e. full infectivity). See https://github.com/jbloomlab/flu_PB1flank-GFP_neut_assay for detailed protocol. 

All neut data relevant to the A/HongKong/45/2019 experiments is documented here. The analysis of A/Perth/16/2009 is limited to visualizing the finalized neut curves and computing the fold-change values. For full analysis of raw neut data from A/Perth/16/2009 assays, see https://github.com/jbloomlab/map_flu_serum_Vietnam_H3_Perth2009/tree/master

* [serum_screening_hk19](serum_screening_hk19): screening for antibodies and sera that neutralize the wildtype library strain HK/19. Includes selection of the 31 sera analyzed in the final H3 HA DMS paper.

* [validations_hk19](validations_hk19): validation of selected mutations by neutralization assays against the wildtype library strain HK/19. Includes MOI tests for each GFP variant virus, analysis of the full set of pilot neutralization assays, and generation of finalized validation neut curves and fold-change-IC50 calculations.

* [validations_perth09](validations_perth09): validation of selected mutations by neutralization assays against the wildtype library strain Perth/09. Includes generation of finalized validation neut curves and fold-change-IC50 calculations.