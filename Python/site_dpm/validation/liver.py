#
# Copyright John Reid 2009
#

"""
Targets of liver specific transcription factors in Wasserman's predictive model for liver paper
"""


from utils import *
from itertools import imap
import biopsy.identifiers.synergizer as S


liver_targets = [
    'G6PC',
    'IGF1',
    'PAH',
    'IGFBP1',
    'CFB', # BF
    'FABP2',
    'GUCA2B',
    'HOXA4',
    'SLC34A1'
]

translated = S.translate('ensembl', 'Mus musculus', 'mgi_symbol', 'ensembl_gene_id', liver_targets)
print translated
print 'Ensembl translated: %d/%d' % (S.how_many_have_translations(translated), len(liver_targets))
ensembl_genes = S.get_translations(translated)
open('liver-targets.txt', 'w').write('\n'.join(ensembl_genes))


liver_ensembl_targets = {
    'G6PC'    : 'ENSMUSG00000078650',
    'IGF1'    : 'ENSMUSG00000020053',
    'PAH'     : 'ENSMUSG00000020051',
    'IGFBP1'  : 'ENSMUSG00000020429',
    'CFB'     : 'ENSMUSG00000024371',
    'FABP2'   : 'ENSMUSG00000023057',
    'GUCA2B'  : 'ENSMUSG00000032978',
    'HOXA4'   : 'ENSMUSG00000000942',
    'SLC34A1' : 'ENSMUSG00000021490',
}





liver_tfs = [
    'HNF1',
    'HNF3',
    'HNF4',
    'cEBP',
]
"""
Liver specific transcription factors as listed in TREMOR paper.
"""

liver_tf_genes = set()
for tf in liver_tfs:
    print tf
    matrices = list(matrices_by_name(tf))
    #print ' '.join(imap(str, ((m, m.name) for m in matrices)))
    genes = set(chain(*(pssm_map[str(m.acc)] for m in matrices)))
    print genes
    print genes.intersection(factors)
    liver_tf_genes.update(genes)
open('liver-tf-ensembl.txt', 'w').write('\n'.join(liver_tf_genes))
