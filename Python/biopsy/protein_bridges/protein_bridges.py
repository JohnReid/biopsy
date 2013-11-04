#
# Copyright John Reid 2007
#

from copy import deepcopy

def points_on_circle( n ):
    from math import pi, sin, cos
    return [
            ( sin( 2 * pi * i / n ), cos( 2 * pi * i / n ) )
            for i in xrange( n )
    ]


def id_to_vertex_map( g, property_name = 'node_id' ):
    "Returns a map from interactors to vertices"
    ids = g.vertex_properties[ property_name ]
    return dict(
            ( ids[ v ], v )
            for v in g.vertices
    )

def get_shortest_distances( g, vertices = None ):
    "The shortest distance between the given vertices and all others in the graph"
    import boost.graph as bgl
    if None == vertices: vertices = g.vertices
    predecessor_map = g.add_vertex_property( 'predecessor', 'vertex' )
    distance_map = g.add_vertex_property( 'distance', 'float' )
    index_map = g.add_vertex_property( 'index', 'index' )
    result = dict( (index_map[ v ], dict()) for v in g.vertices )
    for v1 in vertices:
        if v1 in result: continue # avoid duplicating work
        bgl.dijkstra_shortest_paths(
                g,
                v1,
                predecessor_map,
                distance_map
        )
        for v2 in g.vertices:
            result[ index_map[ v1 ] ][ index_map[ v2 ] ] = distance_map[ v2 ]
    del g.vertex_properties[ 'predecessor' ]
    del g.vertex_properties[ 'distance' ]
    return result

def get_distance_to_vertices( g, v, shortest_distances, vertices ):
    "Get distance from vertex, v, to closest vertex in vertices"
    index_map = g.vertex_properties[ 'index' ]
    return min(
            [
                    shortest_distances[ index_map[ v2 ] ][ index_map[ v ] ]
                    for v2 in vertices
            ]
    )


def _all_paths_recurse(g, u, v, cur_path, paths, max_length=0):
    'Puts all the paths between u and v in the paths argument'
    if max_length and max_length < len(cur_path):
        return
    cur_path.append(u)
    if u == v:
        paths.append(list(cur_path))
    else:
        for w in g.adjacent_vertices(u):
            if w not in cur_path:
                _all_paths_recurse(g, w, v, cur_path, paths, max_length=max_length)
    cur_path.pop()

def all_paths(g, u, v, max_length=0):
    'Returns a list of all the paths between u and v in g'
    paths = []
    _all_paths_recurse(g, u, v, [], paths, max_length=max_length)
    return paths



class protein_bridge_builder( object ):
    """
    Builds a graph of protein-protein interactions that bridge the sets of
    proteins given to it
    """

    def __init__(
            self,
            ppi_network,
            interactor_map
    ):
        """
        Initialises a protein bridge builder.

                - ppi_network is a protein-protein interaction network
                - interactor_map maps remos to interactors in the network
        """
        self.ppi_network = deepcopy(ppi_network)
        self.g = self.ppi_network.g
        self.interactor_map = interactor_map

    def build(
            self,
            verbose = True,
            max_path_length = 1
    ):
        """
        Builds a protein bridge between the remos.
        """
        # Build a map from remos to their vertices
        self.remo_2_vertices_map = dict(
          (
                remo,
                set(self.ppi_network.node_for(i) for i in interactors)
                )
                for remo, interactors
                in self.interactor_map.iteritems()
        )

        # Find the set of all interactors that are associated with at least one remo
        if verbose: print 'Finding proteins associated with remos'
        self.interesting = set()
        for remo, interactors in self.interactor_map.iteritems():
            self.interesting.update(interactors)
        if verbose: print 'Found %d interactors associated with remos' % len( self.interesting )
        if 0 == len( self.interesting ): raise RuntimeError( 'No interactors found for remos' )

        # Calculate the shortest distances between all pairs in the interesting set
        # of interactors
        if verbose: print 'Calculating shortest distances'
        #from IPython.Debugger import Pdb; Pdb().set_trace()
        self.interesting_vertices = [self.ppi_network.node_for(i) for i in self.interesting]
        self.shortest_distances = get_shortest_distances( self.g, self.interesting_vertices )

        # Get rid of those vertices that are too far away from remos to be on a
        # shortest path
        if verbose: print 'Removing proteins too far from remos'
        original_num_vertices = len( self.g.vertices )
        to_remove = [
                v
                for v in self.g.vertices
                if not self.is_inbetween_remos(
                        v,
                        self.shortest_distances,
                        max_path_length
                )
        ]
        if verbose: print 'Removing %d/%d vertices' % ( len( to_remove ), len( self.g.vertices ) )
        for v in to_remove:
            self.g.clear_vertex( v )
            self.g.remove_vertex( v )

        shape = self.g.add_vertex_property( 'shape', 'string' )
        for v in self.g.vertices: shape[ v ] = 'ellipse'

        label = self.ppi_network.add_labels_to_graph()
        url = self.ppi_network.add_urls_to_graph()

        # add vertices for remos
        print 'Adding vertices for remos'
        pos = self.g.add_vertex_property( 'pos', 'string' )
        style = self.g.add_edge_property('style', 'string')
        for e in self.g.edges:
            style[e] = 'solid'
        self.remo_vertices = { }
        for ( (x,y), (remo, vertices) ) in zip(
                points_on_circle( len( self.remo_2_vertices_map ) ),
                self.remo_2_vertices_map.iteritems()
        ):
            v = self.g.add_vertex()
            self.remo_vertices[ remo ] = v
            shape[ v ] = 'diamond'
            label[ v ] = remo
            pos[ v ] = '%f,%f!' % (x * 5,y * 5)
            for v2 in vertices:
                if v2 in self.g.vertices:
                    e = self.g.add_edge( v, v2 )
                    style[e] = 'dotted'

        # remove all leaf nodes (that are not remo nodes)
        while True:
            to_remove = self.leaves_not_in_remos()
            if 0 == len( to_remove ): break
            for v in to_remove:
                self.g.clear_vertex( v )
                self.g.remove_vertex( v )
        if verbose: print 'Left with %d vertices after pruning leaves' % len( self.g.vertices )

        # remove all nodes that are only connected to remos
        to_remove = []
        for v in self.g.vertices:
            for u in self.g.adjacent_vertices(v):
                if u not in self.remo_vertices.values():
                    break
            else:
                to_remove.append(v)
        for v in to_remove:
            self.g.clear_vertex(v)
            self.g.remove_vertex(v)
        if verbose: print 'Left with %d vertices after removing chains of length 1' % len( self.g.vertices )

    def vertices_for_interactors( self, interactors ):
        "Yield all the (interactor, vertex) pairs we have in our id to vertex map"
        for i in interactors:
            if i in self.id_2_vertex:
                yield i, self.id_2_vertex[ i ]

    def leaves_not_in_remos( self ):
        return [
                v
                for v in self.g.vertices
                if 1 == self.g.in_degree( v ) and v not in self.remo_vertices.values()
        ]

    def is_inbetween_remos(
            self,
            v,
            shortest_distances,
            max_path_length = 6
    ):
        import heapq
        heap = []
        # for each remo's vertices find the shortest distance from this vertex
        for remo, vertices in self.remo_2_vertices_map.iteritems():
            if not len(vertices):
                heapq.heappush(heap, max_path_length+1)
            else:
                heapq.heappush(
                        heap,
                        get_distance_to_vertices(
                                self.g,
                                v,
                                shortest_distances,
                                vertices
                        )
                )
        #print heapq.nsmallest( heap, 2 )
        min_distance = sum( heapq.nsmallest( 2, heap ) )
        return min_distance <= max_path_length

    def remo_remo_paths(self, max_length=0):
        'Yield pairs of remos with all paths between each pair'
        if max_length:
            # add 2 to max length to include remo vertices
            max_length += 2
        for remo_1, v1 in self.remo_vertices.iteritems():
            for remo_2, v2 in self.remo_vertices.iteritems():
                yield (remo_1, remo_2), all_paths(self.g, v1, v2, max_length=max_length)

    def path_as_names(self, path):
        ids = self.g.vertex_properties['label']
        return [ ids[v] for v in path[1:-1] ]

    def write_remo_remo_paths(self, filename, max_length=0):
        f = open(filename, 'w')
        for (remo_1, remo_2), paths in self.remo_remo_paths(max_length=max_length):
            for path in paths:
                f.write('%s-%s: %s\n' % (remo_1, remo_2, ' - '.join(self.path_as_names(path))))
            f.write('\n')

    def interactors_for_remo(self, remo):
        "Return the interactors for the given remo that are left in the network."
        return [
          self.ppi_network.interactor_for_node(v)
                for v
                in self.g.adjacent_vertices(self.remo_vertices[remo])
              ]
