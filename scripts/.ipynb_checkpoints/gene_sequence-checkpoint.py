"""Implements ``snakemake`` rule to get protein sequence."""


import Bio.SeqIO

geneseqs = set([])
gene = snakemake.config['gene']
for name, amplicon_file in snakemake.input.items():
    print(f"Reading amplicon from {amplicon_file=}")
    amplicon_of_interest = snakemake.config['amplicons_of_interest'][name]
    amplicon = [a for a in Bio.SeqIO.parse(amplicon_file, 'genbank')
                if a.name == amplicon_of_interest]
    if len(amplicon) != 1:
        raise ValueError(f"{amplicon_file} lacks {amplicon_of_interest=}")
    amplicon = amplicon[0]
    gene_feature = [f for f in amplicon.features if f.type == gene]
    if len(gene_feature) != 1:
        raise ValueError(f"{amplicon_file=} does not have exactly one "
                         f"feature of type {gene=}")
    geneseq = gene_feature[0].extract(amplicon).seq
    print(f"Read amplicon of length {len(geneseq)} from {amplicon_file=}")
    geneseqs.add(geneseq)
if len(geneseqs) != 1:
    raise ValueError(f"Found {len(geneseqs)} different sequences for {gene=}")
geneseq = list(geneseqs)[0]
protseq = geneseq.translate(cds=True)

for outfile, seq in [(snakemake.output.codon, geneseq), (snakemake.output.prot, protseq)]:
    print(f"Writing to {outfile}")
    with open(outfile, 'w') as f:
        f.write(f">{snakemake.config['gene']}\n{str(seq)}\n")
