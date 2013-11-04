#
# Copyright John Reid 2008
#

"""
Code to partition Transfac factors into equivalent sets.
"""

import biopsy.transfac as T
import biopsy
import boost.graph as bgl
from cookbook.lru_cache import lru_cache

class Graph(bgl.Graph):
    __getstate_manages_dict__ = 1
    """So Boost.python knows we manage the object's dict..."""

    def __getstate__(self):
        return (
          bgl.Graph.__getstate__(self),
          self.vertex_id_prop_name,
          self.vertex_id_prop_type
        )

    def __setstate__(self, state):
        bgl_state, self.vertex_id_prop_name, self.vertex_id_prop_type = state
        bgl.Graph.__setstate__(self, bgl_state)
        self.id_2_vertex = dict()
        self.vertex_2_id = self.vertex_properties[self.vertex_id_prop_name]
        for v in self.vertices:
            self.id_2_vertex[self.vertex_2_id[v]] = v

    def __init__(self, vertex_id_prop_name='label', vertex_id_prop_type='string'):
        """
        Creates a new Graph that has a property map from the given type to the vertices.

        @arg vertex_id_prop_name: The name of the property map that maps vertices to ids.
        @arg vertex_id_prop_type: The type of the property map that maps vertices to ids. It can be
        one of the listed types.

           Name         C++ type
        --------------------
           integer      int
           float        float
           vertex       vertex_descriptor
           edge         edge_descriptor
           string       boost::python::str
           point2d      boost::graph::python::point2d
           point3d      boost::graph::python::point3d
           object       boost::python::object
           color        boost::default_color_type
           index        int (contains index of each vertex)
        """
        bgl.Graph.__init__(self)

        self.vertex_id_prop_name = vertex_id_prop_name
        "The name of the property map that maps vertices to their ids."

        self.vertex_id_prop_type = vertex_id_prop_type
        "The type of the property map that maps vertices to their ids (i.e. the type of the ids)."

        self.id_2_vertex = dict()
        "A dict mapping ids to vertices."

        self.vertex_2_id = self.add_vertex_property(vertex_id_prop_name, vertex_id_prop_type)
        "A boost.graph property map mapping vertices to ids."

    def get_id(self, v):
        """
        Return the id for this vertex
        """
        return self.vertex_2_id[v]

    def get_vertex_by_id(self, id):
        """
        Get the vertex with the given id.
        """
        if id not in self.id_2_vertex:
            raise RuntimeError('Id is not in graph.')
        return self.id_2_vertex[id]

    def get_or_add_vertex_by_id(self, id):
        """
        Get the vertex with the given id or if it is not in graph, then add the vertex.
        """
        if id not in self.id_2_vertex:
            v = self.add_vertex()
            self.id_2_vertex[id] = v
            self.vertex_2_id[v] = id
        return self.get_vertex_by_id(id)

    def remove_vertex_by_id(self, id):
        """
        Remove the vertex with the given id from the graph.
        """
        self.remove_vertex(self, self.get_vertex_by_id(id))

    def remove_vertex(self, v):
        """
        Remove the vertex from the graph. Call clear_vertex first if v has edges.
        """
        del self.id_2_vertex[self.vertex_2_id[v]]
        return bgl.Graph.remove_vertex(self, v)

    def add_edge_by_id(self, id1, id2):
        """
        Add an edge between the vertices with the given ids.
        """
        return self.add_edge(self.get_vertex_by_id(id1), self.get_vertex_by_id(id2))

@lru_cache(maxsize=1)
def build_factor_synonyms_graph():
    """
    Build a graph that encodes all the factor synonyms in transfac.
    """
    from itertools import chain
    g = Graph()
    for f in T.Factor.all(): # for each factor
        for synonym1 in chain([f.name], f.synonyms): # for each synonym
            v1 = g.get_or_add_vertex_by_id(synonym1)
            for synonym2 in chain([f.name], f.synonyms): # add an edge to each other synonym
                if synonym1 != synonym2:
                    v2 = g.get_or_add_vertex_by_id(synonym2)
                    if v2 not in g.adjacent_vertices(v1):
                        g.add_edge(v1, v2)
    return g

def remove_small_components(g, num_components, component_map, min_size=2):
    import numpy
    component_sizes = numpy.zeros((num_components,))
    for v in g.vertices:
        component_sizes[component_map[v]] += 1
    for v in g.vertices:
        if component_sizes[component_map[v]] < min_size:
            g.clear_vertex(v)
            g.remove_vertex(v)
    return component_sizes

class FactorSynonyms(object):
    """
    Partitions the set of all factor names into equivalence partitions based on synonyms.

    Maps from factor names to indexes of the partition.
    """

    def __init__(self):
        self.g = build_factor_synonyms_graph()
        self.component_map = self.g.add_vertex_property(name='connected_components', type='integer')
        self.num_components = bgl.connected_components(self.g, self.component_map)
        self._build_partition_synonyms()

    def _build_partition_synonyms(self):
        """
        Calculates one synonym to represent each partition
        """
        self.partition_synonyms = [None] * self.num_components
        for v in self.g.vertices:
            idx = self.component_map[v]
            if None == self.partition_synonyms[idx]:
                self.partition_synonyms[idx] = self.g.get_id(v)

    def get_partition_idx(self, factor_name):
        """
        Get the index of the partition that this factor name is in.
        """
        v = self.g.get_vertex_by_id(factor_name)
        return self.component_map[v]

    def get_partition_synonym(self, partition_idx):
        """
        Return the representative synonym for this partition
        """
        return self.partition_synonyms[partition_idx]

    def get_partition_synonyms(self, partition_idx):
        """
        Return the synonyms that make up this partition
        """
        return [
          self.g.get_id(v)
          for v in self.g.vertices
          if partition_idx == self.component_map[v]
        ]

    def get_synonym(self, factor_name):
        """
        Return the representative synonym of this factor name
        """
        return self.get_partition_synonym(self.get_partition_idx(factor_name))

    def get_synonyms(self, factor_name):
        """
        Return all the synonyms of this factor name
        """
        return self.get_partition_synonyms(self.get_partition_idx(factor_name))

class Pssm2FactorSynonymMap(dict):
    """
    Maps Transfac PSSM accessions to sets of factor synonyms
    """

    def __init__(self, factor_synonyms):
        self.factor_synonyms = factor_synonyms
        for acc in biopsy.get_transfac_pssm_accessions(biopsy.transfac.PssmFilter.all_pssms()):
            for factor in biopsy.transfac.TableLink(acc).entry.factors:
                self[acc].add(self.factor_synonyms.get_synonym(factor.link.entry.name))

    def __missing__(self, k):
        self[k] = set()
        return self[k]

if '__main__' == __name__:
    factor_synonyms = FactorSynonyms()
