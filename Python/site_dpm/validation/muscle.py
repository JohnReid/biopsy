#
# Copyright John Reid 2009
#


muscle_tfs = [
    'MYOD',
    'SRF',
    'MEF2',
    'TEF1',
    'SP1',
]
"""
Muscle specific transcription factors as listed in TREMOR paper.
"""

muscle_tf_genes = set()
for tf in muscle_tfs:
    print tf
    matrices = list(matrices_by_name(tf))
    #print ' '.join(imap(str, ((m, m.name) for m in matrices)))
    genes = set(chain(*(pssm_map[str(m.acc)] for m in matrices)))
    print genes
    print genes.intersection(factors)
    muscle_tf_genes.update(genes)
open('muscle-tf-ensembl.txt', 'w').write('\n'.join(muscle_tf_genes))
