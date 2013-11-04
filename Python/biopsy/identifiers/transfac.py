#
# Copyright John Reid 2007
#

"""
Code to help map transfac identifiers onto entrez genes/proteins, swissprot ids
"""

import biopsy, csv, cookbook, itertools

def is_mouse_gene(g):
    for r in g.db_refs:
        if 'ENSMUSG' == r.table or biopsy.transfac.db.mgi == r.db:
            return True
    return False

def to_mouse_map(ortho_filename = os.path.join(biopsy.get_data_dir(), 'TreeFam', 'orthologs', 'MOUSE_ORTHO.tsv')):
    return dict(
            (biopsy.DbRef.try_to_parse(row[1]), biopsy.DbRef.try_to_parse(row[0]))
            for row
            in csv.reader(
                    open(ortho_filename),
                    delimiter = '\t'
            )
            if len(row) > 1
    )

def refs_for_factor(f):
    for r in f.db_refs:
        yield r
    if f.gene:
        for r in f.gene.entry.db_refs:
            yield r

def refs_for_pssm(m):
    for f_ref in m.factors:
        f = f_ref.link.entry
        if f_ref.species.startswith('mouse') or (f.gene and is_mouse_gene(f.gene.entry)):
            for r in refs_for_factor(f):
                yield r

def names_for_pssm(m):
    for f_ref in m.factors:
        f = f_ref.link.entry
        yield f.name
        for s in f.synonyms:
            yield s

if '__main__' == __name__:
    import biopsy.transfac as T

    matrix_references = cookbook.DictOfLists()
    matrix_names = cookbook.DictOfLists()

    for m in T.Matrix.all():
        acc = m.acc
        for r in refs_for_pssm(m):
            matrix_references[acc].append(r)
        for n in names_for_pssm(m):
            matrix_names[acc].append(n)

    #for m in T.Site.all():
    #       refs = [r for r in refs_for_pssm(m)]
    #       if len(refs):
    #               references[m] = refs

    print 'Found references for %4d of %4d matrices' % (len(matrix_references), len(T.Matrix.all()))
