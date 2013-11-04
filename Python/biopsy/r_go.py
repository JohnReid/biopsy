#
# Copyright John Reid 2006
#

import rpy

rpy.r.library("Category")
rpy.r.library("biomaRt")

rpy.r.source('../../R/go_categorise.R')

_marts = { }

def get_mart( dataset = "mmusculus_gene_ensembl" ):
    """Gets ensembl mart for dataset"""
    if not dataset in _marts:
        _marts[ dataset ] = rpy.with_mode(
                rpy.NO_CONVERSION,
                rpy.r.useMart
        )(
                "ensembl",
                dataset = dataset
        )
    return _marts[ dataset ]

def categorise_genes( mart, genes ):
    """Annotate/categorise genes according to GO"""
    import rpy
    return rpy.with_mode(
            rpy.NO_CONVERSION,
            rpy.r.categoriseGenes
    )(
            genes,
            mart
    )
