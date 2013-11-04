#
# Copyright John Reid 2007
#

import biopsy, os, csv, itertools
import biopsy.transfac as trans

def one_to_many( relationships ):
    """Takes a sequence of binary relationships and creates a dict
    mapping from one-to-many"""
    d = { }
    for o1, o2 in relationships:
        if d.has_key( o1 ):
            d[o1].add( o2 )
        else:
            d[o1] = set( ( o2, ) )
    return d

def collapse_relationships( relationships ):
    """Takes a sequence of dicts representing one-to-many relationships and
    collapses this to one dict representing overall one-to-many relationships"""
    overall = { }
    for k, v in relationships[0]:
        pass

def take_first_2( rows ):
    """Yield first 2 elements in each row as a tuple"""
    for row in rows:
        try: yield ( row[0], row[1] )
        except: pass

def convert_to_db_refs( rows ):
    """Convert each element in a row to an DbRef"""
    for row in rows:
        try: yield [ biopsy.DbRef( acc ) for acc in row ]
        except: pass


def generate_db_refs_for_mouse_ortho( reader ):
    for row in convert_to_db_refs( take_first_2( reader ) ):
        yield (row[1], row[0])

# load mouse ortholog data
def load_mouse_ortho_data(
        filename = os.path.join(biopsy.get_data_dir(), 'TreeFam', 'orthologs', 'MOUSE_ORTHO.tsv')
):
    reader = csv.reader(
            open( filename, "rb" ),
            delimiter = '\t'
    )
    return one_to_many( generate_db_refs_for_mouse_ortho( reader ) )

_mouse_ortho_data = None
def mouse_ortho_data():
    global _mouse_ortho_data
    if None == _mouse_ortho_data:
        _mouse_ortho_data = load_mouse_ortho_data()
    return _mouse_ortho_data

def print_mouse_ortho_stats( ortho ):
    # print mouse ortho data stats
    print 'Mapped %d identifiers to mouse genes' % len( ortho )
    multiple_genes = 0
    mouse_universe = set()
    for g, s in ortho.iteritems():
        if len( s ) > 1: multiple_genes += 1
        mouse_universe.update( s )
    print 'Mapped %d identifiers to multiple mouse genes' % multiple_genes
    print 'Mapped onto %d mouse genes' % len( mouse_universe )

def ensembl_mouse_for_ensembl_ref( ref ):
    "Yields ensembl mouse genes for given ensembl reference"
    assert trans.db.ensembl == ref.db
    if 'ENSMUSG' == ref.table: yield ref # it is already a mouse gene
    if mouse_ortho_data().has_key( ref ):
        for g in mouse_ortho_data()[ ref ]:
            yield g

def ensembl_mouse_for_ref( ref ):
    "Yields ensembl mouse genes for reference"
    # is it an ensembl reference?
    if trans.db.ensembl == ref.db:
        for g in ensembl_mouse_for_ensembl_ref( ref ):
            yield g

def ensembl_mouse_for_transfac_pssm( p ):
    "Yields ensembl mouse genes for transfac pssm"
    for f in p.factors:
        g = f.link.entry.gene
        if g.known:
            for r in g.entry.db_refs:
                for em in ensembl_mouse_for_ref( r ):
                    yield em

def encoding_genes_for_transfac_pssms( pssms ):
    for p in pssms:
        yield (p.acc, [ g for g in ensembl_mouse_for_transfac_pssm( p ) ])



if "__main__" == __name__:
    print_mouse_ortho_stats( mouse_ortho_data() )
    genes_for_transfac_pssms = dict(
            encoding_genes_for_transfac_pssms(
                    itertools.chain(
                            trans.Site.pssms(),
                            trans.Matrix.pssms()
                    )
            )
    )
