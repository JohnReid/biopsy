#
# Copyright John Reid 2009
#

"""
Code to generate validation sets from the literature, KEGG etc..
"""



from shared import *
import logging


def get_literature_sets():
    logging.info('Getting literature validation sets.')
    import literature
    return literature.map_factor_sets(literature.tf_sets_by_name)

@global_cached_method('kegg-sets')
def get_kegg_sets():
    logging.info('Getting KEGG validation sets.')
    import kegg as K
    return dict(
        ('kegg-%s-%s' % (id, K.get_kegg_pathway_info(id).split(';')[0]), targets)
        for id, targets in K.get_kegg_pathways()
    )


def get_symatlas_sets():
    logging.info('Getting SymAtlas validation sets.')
    import biopsy.data.symatlas as SA
    probes_to_genes = SA.probes_to_genes()
    dataset, tissues, probes, fold_changes = SA.expression_data()
    highly_expressed = SA.fold_change_above_median(fold_changes, fold_change=options.symatlas_fold_change_threshold)
    highly_expressed_probes = SA.match_probe_sets(highly_expressed.T, probes)
    highly_expressed_genes = [
        set(probes_to_genes[p] for p in hep if p in probes_to_genes and probes_to_genes[p])
        for hep in highly_expressed_probes
    ]
    return dict(zip(('symatlas-%s' % t for t in tissues), highly_expressed_genes))


@log_exceptions()
@caching_decorator('validation_sets')
def generate_validation_sets():
    """
    Return 2 dicts mapping validation set names to their ensembl IDs. One for factors, one for targets.
    """
    factor_sets = dict()
    target_sets = dict()

    kegg_sets = get_kegg_sets()
    logging.info('Have %d validation sets from KEGG.', len(kegg_sets))
    target_sets.update(kegg_sets)

    symatlas_sets = get_symatlas_sets()
    logging.info('Have %d validation sets from SymAtlas.', len(symatlas_sets))
    target_sets.update(symatlas_sets)

    literature_sets = get_literature_sets()
    logging.info('Have %d factor validation sets from literature.', len(literature_sets))
    factor_sets.update(literature_sets)

    return factor_sets, target_sets



if '__main__' == __name__:
    generate_validation_sets()
