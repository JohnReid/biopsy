#
# Copyright John Reid 2009
#

"""
Code to use mouse promoters from TransPro
"""

import biopsy.data.biobase.transpro as transpro
import biopsy.identifiers.biomart as biomart
from analysis import *
from cookbook.dicts import DictOf


def has_ref_fn(ref_type):
    def has_ref(promoter):
        for r in promoter.refs:
            if r.startswith(ref_type):
                return True
        return False
    return has_ref

def mgi_id_for(promoter):
    for r in promoter.refs:
        if r.startswith('MGI'):
            return r.split(':')[1].split(';')[0].strip()
    return None

def get_mgi_to_ensembl_map():
    return dict(
      biomart.quick_query(
        dataset='mmusculus_gene_ensembl',
        attributes=['mgi_id', 'ensembl_gene_id'],
      )
    )

def ensembl_for(promoter):
    return mgi_to_ensembl.get('MGI:%s' % mgi_id_for(promoter), None)


@global_cached_method('entrez-to-ensembl')
def entrez_to_ensembl():
    import biopsy.identifiers.biomart as B
    return dict(
        B.quick_query(
            dataset='mmusculus_gene_ensembl',
            attributes=['entrezgene', 'ensembl_gene_id'],
            filters=()
        )
    )

if '__main__' == __name__:
    mouse_promoters = transpro.get_mouse_promoters()

    for ref_type in [
        'MGI',
        'ENTREZGENE',
        'REFSEQ',
        'UNIGENE',
        'ENSMUSG'
    ]:
        print ref_type, len(filter(has_ref_fn(ref_type), mouse_promoters))

    print map(mgi_id_for, mouse_promoters[:20])

    mgi_to_ensembl = get_mgi_to_ensembl_map()

    ensembl_promoters = DictOf(list)
    for p in mouse_promoters:
        ensembl = ensembl_for(p)
        if ensembl:
            ensembl_promoters[ensembl].append(p)

    sequence_analyser = get_sequence_analyser()
    analysis = DictOf(list)
    for ensembl, remos in ensembl_promoters.iteritems():
        for remo in remos:
            analysis[ensembl].append(sequence_analyser(remo.sequence))


    fasta = open('mouse-promoters.fa', 'w')
    for ensembl, remos in ensembl_promoters.iteritems():
        for i, remo in enumerate(remos):
            fasta.write('> %s - %d\n' % (ensembl, i))
            fasta.write(remo.sequence)
            fasta.write('\n')
    fasta.close()
