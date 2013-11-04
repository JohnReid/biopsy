#
# Copyright John Reid 2008
#


"""
Code to access KEGG data and test some gene sets for enrichment.
"""

import biopsy, os
from cookbook.cache_decorator import cachedmethod


@cachedmethod
def get_wsdl_proxy(wsdl):
    from SOAPpy import WSDL
    return WSDL.Proxy(wsdl)

@cachedmethod
def kegg_service():
    return get_wsdl_proxy('http://soap.genome.jp/KEGG.wsdl')

@cachedmethod
def get_kegg_pathway_info(pathway):
    return kegg_service().btit('path:%s' % pathway)[14:].strip()

def get_kegg_pathways():
    pathway_dir = os.path.join(biopsy.get_data_dir(), 'KEGG', 'pathways')
    return [
      (file, set(l.strip() for l in open(os.path.join(pathway_dir, file))))
      for file in os.listdir(pathway_dir)
    ]


def as_set(s):
    if set == type(s):
        return s
    else:
        return set(s)

def as_str_set(s):
    return as_set(str(t) for t in s)


if '__main__' == __name__:
    from gene_set_enrichment import test_enrichment

    kegg_pathways = get_kegg_pathways()

    for i, tp in enumerate(transcriptional_programs):
        threshold = 0.05
        target_universe = as_str_set(tp.genes)
        factor_universe = as_str_set(tp.factors)
        target_set = as_str_set(tp.tp_targets)
        factor_set = as_str_set(tp.tp_factors)
        for name, pathway in kegg_pathways:
            pathway_target_set = as_str_set(pathway)
            pathway_target_set.intersection_update(target_universe)
            pathway_factor_set = as_str_set(pathway)
            pathway_factor_set.intersection_update(factor_universe)
            targets_p_value = test_enrichment(target_set, pathway_target_set, target_universe) * len(kegg_pathways)
            factors_p_value = test_enrichment(factor_set, pathway_factor_set, factor_universe) * len(kegg_pathways)
            if targets_p_value < threshold or factors_p_value < threshold:
                print i, name, targets_p_value, factors_p_value
