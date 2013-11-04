#
# Copyright John Reid 2009
#

"""
Code to examine factors in Transfac.
"""

import biopsy.transfac as T

def factors_by_name(factor_name):
    "@return: Yields factors that have the given name or are synonymous with it."
    factor_name = factor_name.upper()
    for f in T.Factor.all():
        if f.name.upper() == factor_name:
            yield f
            continue
        for synonym in f.synonyms:
            if synonym.upper() == factor_name:
                yield f
                break

def matrices_for_factor(factor):
    "@return: Yields matrices for the factor."
    for p in T.Matrix.all():
        for f in p.factors:
            if factor.acc == f.link:
                yield p
                break

def matrices_by_factor_name(factor_name):
    "@return: Yields matrices that can be linked to the factor_name using factors_by_name and matrices_for_factor."
    for f in factors_by_name(factor_name):
        for p in matrices_for_factor(f):
            yield p


def ens_genes_for(factor_name, pssm_map):
    """
    The ensembl genes for the factor name (using the given pssm map).
    """
    return set(
        pssm_map[m]
        for m
        in set(str(m.acc) for m in matrices_by_factor_name(factor_name))
        if m in pssm_map
    )


def map_factor_names_to_ensembl(factor_names):
    "@return: A map from the factor names to ensembl genes using the given map."
    return dict((tf, ens_genes_for(tf, pssm_map)) for tf in factor_names)





if '__main__' == __name__:
    e2f_factors = list(factors_by_name('e2f'))
    clox_factors = list(factors_by_name('clox'))

    e2f_matrices = set(matrices_by_factor_name('e2f'))
    clox_matrices = set(matrices_by_factor_name('clox'))
