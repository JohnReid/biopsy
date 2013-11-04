#
# Copyright John Reid 2006
#

import biopsy.transfac as trans

import itertools

def print_sizes():
    print '# factors:   %d' % len( trans.Factor.all() )
    print '# genes:     %d' % len( trans.Gene.all() )
    print '# sites:     %d' % len( trans.Site.all() )
    print '# fragments: %d' % len( trans.Fragment.all() )
    print '# matrices:  %d' % len( trans.Matrix.all() )
# print_sizes(); raise

def factors_for( p ):
    "Yield all the factors for the pssm"
    for f in p.factors:
#               try:
        yield f.link.entry
#               except: pass

def genes_for( p ):
    "Yield all the genes for the pssm that have an ensembl reference"
    for f in factors_for( p ):
        g = f.gene
        if g.known and gene_has_ensembl_ref( g.entry ):
            yield g

def gene_has_ensembl_ref( g ):
    "Does the gene have an ensembl reference?"
    for r in g.db_refs:
        if trans.db.ensembl == r.db: return True
    return False

genes_for_pssms = dict(
        [
                ( p.acc, [ g for g in genes_for( p ) ] )
                for p in itertools.chain( trans.Matrix.pssms(), trans.Site.pssms() )
        ]
)
print '# pssms:', len( genes_for_pssms )
print '# pssms with ensembl genes:', len(
        [
                a
                for a, g in genes_for_pssms.iteritems()
                if len(g)
        ]
)
