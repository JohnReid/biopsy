#
# Copyright John Reid 2007
#

"""
Code to build a map between various different identifier types
"""

# Import access to various databases
from . import biomart
from . import ensembl
from . import entrez
from . import transfac
from . import uniprot

from cookbook.persisted_cache import PersistedCache

import biopsy
T = biopsy.transfac

class IdentifierMap(object):
    """
    A map between various database identifiers.

    Stores relationships in a boost.graph.
    """

    def __init__(self):
        import boost.graph as bgl
        self.g = bgl.Digraph()
        self.ids = self.g.add_vertex_property('id', type='object') # map from vertices to ids
        self.vertices = dict() # dict from id to vertices

    def __getinitargs__():
        return tuple()

    def __getstate__(self):
        vertices = set(self.ids[v] for v in self.g.vertices)
        edges = set((self.ids[self.g.source(e)], self.ids[self.g.target(e)]) for e in self.g.edges)
        return (vertices, edges)

    def __setstate__(self, t):
        self.__init__()
        vertices, edges = t
        for v in vertices:
            self.get_vertex(v)
        for source, target in edges:
            self.add_edge(source, target)

    def vertex(self, id):
        "Returns the vertex for the id, does not create new if not there"
        if id not in self.vertices:
            return None
        return self.vertices[id]

    def get_vertex(self, id):
        "Returns the vertex for the id, creates new if not there"
        if id not in self.vertices:
            v = self.g.add_vertex()
            self.vertices[id] = v
            self.ids[v] = id
        return self.vertices[id]

    def links(self, id, predicate=None):
        """
        Return all the ids connected to the given id that satisfy the predicate
        if given.
        """
        import boost.graph as bgl
        v = self.vertex(id) # is the id in the graph
        vertices = []
        if None != v:
            class PredicateAppendingVisitor(bgl.bfs_visitor):
                def examine_vertex(visitor, v, g):
                    id = self.ids[v]
                    if None == predicate or predicate(id):
                        vertices.append(id)
            bgl.breadth_first_search(self.g, v, visitor=PredicateAppendingVisitor())
        return vertices

    def add_edge(self, id1, id2):
        """
        Add an edge between the ids if not there already. Returns edge.
        """
        v1, v2 = self.get_vertex(id1), self.get_vertex(id2)
        return self.g.add_edge(v1, v2)


def transfac_mapper(r):
    """
    Takes a transfac database reference and yields related entries
    """
    link = T.TableLink.from_db_ref(r)

    # Matrix
    if T.trans_data.matrix == link.db:
        for f_ref in link.entry.factors:
            yield f_ref.link.as_db_ref()

    # Factor
    elif T.trans_data.factor == link.db:
        for ref in link.entry.db_refs:
            yield ref
        if link.entry.gene:
            yield link.entry.gene.as_db_ref()
        for sub_family in link.entry.sub_families:
            yield sub_family.as_db_ref()
        for sub_unit in link.entry.subunits:
            yield sub_unit.as_db_ref()

    # Gene
    elif T.trans_data.gene == link.db:
        for ref in link.entry.db_refs:
            yield ref

    # Site
    elif T.trans_data.site == link.db:
        for f_ref in link.entry.factors:
            yield f_ref.link.as_db_ref()


def ensembl_mapper(r):
    # get mouse orthologs if in suitable other species
    if r.table != 'ENSMUSG':
        if (
          'ENSBTAG' == r.table
          or 'ENSCAFG' == r.table
          or 'ENSDARG' == r.table
          or 'ENSG' == r.table
          or 'ENSGALG' == r.table
          or 'ENSMODG' == r.table
          or 'ENSMUSG' == r.table
          or 'ENSRNOG' == r.table
          or 'ENSXETG' == r.table
          or 'SINFRUG' == r.table
        ):
            for mouse_ortholog in ensembl.orthologs()[(r.table, 'ENSMUSG')][r]:
                yield mouse_ortholog
    else: # yield the proteins associated with the mouse gene
        if r in ensembl.mgi_ids():
            yield ensembl.mgi_ids()[r]
        for transcript, protein in ensembl.proteins()['ENSMUSG'][r]:
            yield protein

def uniprot_mapper(r):
    # get ensembl references
    if r in uniprot.acc2ref():
        for ref in uniprot.acc2ref()[r]:
            yield ref

def entrez_gene_mapper(r):
    "Maps entrez genes to proteins"
    gene_id = r.acc
    if gene_id in entrez.mouse_proteins().geneid2protein:
        for protein_id in entrez.mouse_proteins().geneid2protein[gene_id]:
            yield T.DbRef(T.db.entrez_protein, "", protein_id)

class IdentifierMapper(dict):
    def __init__(self, map):
        self.map = map
    def __missing__(self, db):
        return lambda x: []
    def __call__(self, r):
        self.map.get_vertex(r) # make sure we have our own vertex
        for related in self[r.db](r):
            if r == related:
                continue
            #print r, related
            if None == self.map.vertex(related): # if the map hasn't seen the related ref
                self(related) # add it
            self.map.add_edge(r, related)


def default_map_and_mapper():
    map = IdentifierMap()

    mapper = IdentifierMapper(map)
    mapper[T.db.transfac] = transfac_mapper
    mapper[T.db.ensembl] = ensembl_mapper
    mapper[T.db.entrez_gene] = entrez_gene_mapper
    mapper[T.db.swissprot] = uniprot_mapper

    return map, mapper

def build_map():
    map, mapper = default_map_and_mapper()
    for m in T.Matrix.all():
        mapper(m.acc.as_db_ref())
    for s in T.Site.all():
        if s.id.factor != 'CONS':
            continue
        mapper(s.acc.as_db_ref())
    return map

identifier_map = PersistedCache(build_map, os.path.join(biopsy.get_data_dir(), 'identifiers', 'identifier_map.pickle'))

def small_test_map():
    'Return a small map for testing purposes.'
    map, mapper = default_map_and_mapper()
    mapper(T.DbRef.parse_as('71431', T.db.entrez_gene))
    return map


def matrices_that_map_to(map, db):
    """
    Return a set of those matrices that map to the at least one entry in the
    given database type
    """
    matrices = set()
    for m in T.Matrix.all():
        links = map.links(m.acc.as_db_ref())
        for l in links:
            if l.db == db:
                matrices.add(m)
                break
    return matrices
