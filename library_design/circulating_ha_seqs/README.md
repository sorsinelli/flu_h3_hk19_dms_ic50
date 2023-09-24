# Recently circulating H3 strains

The barcoded H3 library also allows us to spike in full-length strains of interest and analyze neutralization of these strains by different serum samples. We initially identified 13 strains, not including the library backbone A/HongKong/45/2019, that represent currently circulating clades of interest. We then incorporated an additional set of 8 strains that the Hensley lab is analyzing for neutralization.

![H3 Strain Summary](recently_circulating_h3.pdf)

`gen_barcodes.ipynb` was used to generate random barcodes to attach to each strain for synthesis. Strains were then ordered from Twist as dsDNA fragments with overlaps in the first 19 amino acids from WSN, and the read1 sequence downstream of the gene, for direct ligation into the chimeric backbone. 