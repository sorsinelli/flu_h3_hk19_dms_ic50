# Input data

## PacBio full-length variant sequencing
[pacbio_runs.csv](pacbio_runs.csv) contains information on PacBio runs linking the barcodes to the mutations. It has the following columns:
 - `library`: name of the library sequenced
 - `run`: date of the pacbio library submission (use this date to refer to experimental notebook *TODO: link notebooks*)
 - `subreads`: BAM file from running CCS

The amplicon itself is defined in the following files:
 - [hk19_pacbio_amplicon.gb](hk19_pacbio_amplicon.gb): the amplicon for the digested plasmid.
 - [feature_parse_specs.yaml](feature_parse_specs.yaml): how to parse features from the plasmid amplicon with `alignparse`.
 
[site_map.csv](site_map.csv) contains site numbering conversion, and defines the wildtype amino acid sequence for the gene of interest.
 - 'sequential_site': sequential chimeric site numbering, where '1' is the first amino acid in the chimeric construct.
 - 'reference_site': standard H3 numbering, where '1' is codon 20 in the chimeric construct.

## Illumina barcode sequencing
[barcode_runs.csv](barcode_runs.csv) has the Illumina barcode-sequencing runs used to count barcodes in different conditions.
It describes the samples, which should be named as clearly as possible.
It has the following columns:

 - *date*: date experiment was performed in `YYYY-MM-DD` encoding.
 - *virus_batch*: batch of virus used for the experiment.
 - *library*: which virus library was used.
 - *sample_type*: can be one of the following:
   + *no-antibody_control*: entry mediated by HA
   + *antibody*: encompasses sera and antibodies
   + *other_control*: some control unrelated to antibody selections.
 - *antibody*: name of the antibody if this sample has *sample_type* of *antibody*
 - *antibody_concentration*: concentration of antibody if this sample has *sample_type* of antibody. For sera, should be a fraction < 1 giving dilution (**not** a dilution factor).
 - *replicate*: experimental replicate.
 - *fastq_R1*: path to R1 FASTQ file, or semi-colon de-limited list of multiple FASTQs
 - *notes*: any other notes about the sample.

The file [neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) has the barcodes on the neutralization standard that is not expected to be neutralized by the sera or antibodies.
