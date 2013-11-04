#
# Copyright John Reid 2006
#

import biopsy.transfac as trans


def check_for_unknown_db_refs(entries):
    "Checks entries for unknown db_refs"
    num_unknown = 0
    for s in entries:
        for r in s.db_refs:
            if trans.db.unknown_db == r.db:
                num_unknown += 1
    return num_unknown


def print_unknown_db_refs():
    "Prints summary of how many unknown db_refs there are in transfac"
    print '# sites with unknown db_refs: %d' % check_for_unknown_db_refs(trans.Site.all())
    # print '# matrices with unknown db_refs: %d' % check_for_unknown_db_refs(trans.Matrix.all())
    print '# factors with unknown db_refs: %d' % check_for_unknown_db_refs(trans.Factor.all())
    print '# genes with unknown db_refs: %d' % check_for_unknown_db_refs(trans.Gene.all())


def generate_pssm_to_factor_relationships(pssms):
    "Generate pairs of relationships between pssms and factors"
    for p in pssms:
        for f in p.factors:
            if f.link.known:
                yield p.acc, f.link


def pssm_to_factor_map(pssms):
    return one_to_many(generate_pssm_to_factor_relationships(pssms))


def generate_factor_to_gene_relationships(factors):
    for f in factors:
        g = f.gene
        if g.known:
            yield f.acc, g


site_to_factor = pssm_to_factor_map(trans.Site.pssms())
matrix_to_factor = pssm_to_factor_map(trans.Matrix.pssms())
factor_to_gene = one_to_many(generate_factor_to_gene_relationships(trans.Factor.all()))

raise

