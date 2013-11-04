#
# Copyright John Reid 2009
#

import biopsy.identifiers.synergizer as S
import biopsy.identifiers.biomart as B
import csv, os, sys
from itertools import chain, imap
import transfac_factors


def ens_genes_for(factor_name, pssm_map):
    return set(
        pssm_map[m]
        for m
        in set(str(m.acc) for m in transfac_factors.matrices_by_factor_name(factor_name))
        if m in pssm_map
    )

#pssm_map = get_pssm_to_ensembl_map_min_range()

elkon_tfs_small = [
    'E2F',
    'NF-Y',
    'NRF-1',
    'CREB',
]
"""
The transcription factors in "Genome-Wide In Silico Identification of Transcriptional
Regulators Controlling the Cell Cycle in Human Cells" Elkon et al. Table 1
"""




elkon_tfs = [
    'E2F',
    'Sp1',
    'NF-Y',
    'NRF-1',
    'ETF',
    'ATF',
    'CREB',
    'Arnt',
    'YY1',
]
"""
The transcription factors in "Genome-Wide In Silico Identification of Transcriptional
Regulators Controlling the Cell Cycle in Human Cells" Elkon et al. Table 3
"""
elkon_tf_ens_genes = dict((tf, ens_genes_for(tf, pssm_map)) for tf in elkon_tfs)
raise
for factor_name in elkon_tfs:
    print factor_name, ':', ','.join(set(m.name for m in transfac_factors.matrices_by_factor_name(factor_name)))




open('elkon-cell-cycle-small-ensembl.txt', 'w').write('\n'.join('\n'.join(pssm_map[matrix]) for matrix in elkon_tfs_small.values()))
open('elkon-cell-cycle-ensembl.txt', 'w').write('\n'.join('\n'.join(pssm_map[matrix]) for matrix in elkon_tfs.values()))






#translated = S.translate('ensembl', 'Homo sapiens', 'unigene', 'ensembl_gene_id', [ 'Hs.239', 'Hs.715518', 'Hs.181768', 'Hs.427236' ])
#print 'Translated: %s' % ', '.join('->'.join(map(str, x)) for x in translated)

def yield_refs(filename):
    """Yield references in filename"""
    f = open(file)
    f.next() # ignore header line
    for row in csv.reader(f, delimiter='\t'):
        # print row
        yield row[0]
    f.close()

filenames = [
    '1134_and_Cho.txt',
    '1134_not_Cho.txt',
    'cho_not_1134.txt',
]

references = dict()
for file in filenames:
    print '****************** %s ********************' % file
    references[os.path.splitext(file)[0]] = set(yield_refs(file))


# get all the HS unigenes from all the files
cell_cycle_hs_unigene = list(chain(*references.values()))


# translate them to ensembl HS genes using the Synergizer
translated = S.translate('ensembl', 'Homo sapiens', 'unigene', 'ensembl_gene_id', cell_cycle_hs_unigene)
print 'Ensembl translated: %d' % S.how_many_have_translations(translated)
ensembl_hs_genes = S.get_translations(translated)

def yield_mouse_orthologs(hs_genes):
    # map into mouse orthologs using biomart
    query = B.new_query()
    dataset = B.add_dataset(query, 'hsapiens_gene_ensembl')
    B.add_attribute(dataset, 'ensembl_gene_id')
    B.add_attribute(dataset, 'mouse_ensembl_gene')
    filter = B.add_filter(dataset, name='ensembl_gene_id', value='')
    filter.set('value', ','.join(ensembl_hs_genes))
    for chunk in B.split_big_list(ensembl_hs_genes, 50):
        #logging.info('Querying Ensembl biomart for chunk of %d genes', len(chunk))
        filter.set('value', ','.join(chunk))
        for row in B.yield_csv_query_results(query):
            if row[1]:
                yield row[1]
ensembl_mm_genes = set(yield_mouse_orthologs(ensembl_hs_genes))
open('mouse_cell_cycle.txt', 'w').write('\n'.join(ensembl_mm_genes))
