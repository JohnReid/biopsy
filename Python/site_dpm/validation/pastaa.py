#
# Copyright John Reid 2009
#

"""
Sets of TFs defined in PASTAA paper (Vingron)
"""

import biopsy.transfac as T
from cookbook.dicts import DictOf

def get_matrix_by_name(name):
    for m in T.Matrix.all():
        if -1 != m.name.lower().find(name.lower()):
            return m
    return None

pastaa_matrices = {
    'muscle'    : [ 'SRF_Q5_01', 'SRF_01', 'SRF_Q5_02', 'SRF_C', 'MTATA_B' ],
    'heart'     : [ 'MEF2_Q6_01', 'SRF_C', 'RSRFC4_01', 'MTATA_B', 'MEF2_02' ],
    'liver'     : [ 'HNF4_Q6_01', 'HNF1_01', 'HNF4_01', 'HNF1_Q6', 'HNF1_C', ],
    'retina'    : [ 'CRX_Q4', 'CHX10_01' ],
    'leukocyte' : [ 'NFKAPPAB65_01', 'NFKAPPAB_01', 'NFKB_Q6_01', 'CREL_01', 'ETS_Q6' ]
}

pssm_map = get_pssm_to_ensembl_map_min_range()
pastaa_gene_sets = DictOf(set)

for group, names in pastaa_matrices.iteritems():
    for name in names:
        m = get_matrix_by_name(name)
        if m:
            ensembl = pssm_map.get(str(m.acc), '<unknown>')
            if str(m.acc) in pssm_map:
                pastaa_gene_sets[group].add(ensembl)
            logging.info(
                '%20s : %16s = %s (%16s) - %s (%s)',
                group,
                name,
                m.acc,
                m.name,
                ensembl,
                ensembl_names.get(ensembl, '<unknown>')
            )
        else:
            logging.info('%20s : %10s : NO MATCH', group, name)
