#
# Copyright John Reid 2009
#

"""
Gets promomters from results that have a certain GO annotation.
"""

import csv, sys, go, ucsc_promoters, corebio.seq, corebio.seq_io.fasta_io

def get_genes_for_go_in(go_category, program_targets_file):
    for row in csv.reader(open(program_targets_file), delimiter=','):
        if go_category in go_annotations[row[0]]:
            yield row[0]

if '__main__' == __name__:
    go_annotations = go.get_all_ensembl_go_annotations()
    go_category = 'GO:0004984'
    program_targets_file = sys.argv[1]
    genes = set(get_genes_for_go_in(go_category, program_targets_file))
    promoters = ucsc_promoters.get_promoters()
    fasta = open('genes-%s.fa' % go_category, 'w')
    for gene in genes:
        for s, seq in enumerate(promoters[gene]):
            seq = corebio.seq.dna(seq)
            seq.name = '%s - %d' % (gene, s)
            corebio.seq_io.fasta_io.writeseq(fasta, seq)
    fasta.close()
