#
# Copyright John Reid 2007
#


def test_pickle():
    import boost.graph as bgl
    g = bgl.Graph()
    p = g.add_vertex_property( 'prop', 'float' )
    pickle.dump( g, open( 'graph.pickle', 'w' ) )
    pickle.load( open( 'graph.pickle', 'r' ) )


def test_dijkstra():
    import boost.graph as bgl
    graph = bgl.Graph()
    a = graph.add_vertex()
    b = graph.add_vertex()
    e = graph.add_edge(a, b)
    distances = graph.add_vertex_property('float')
    predecessors = graph.add_vertex_property('vertex')
    bgl.dijkstra_shortest_paths(
            graph,
            a,
            predecessor_map = predecessors,
            distance_map = distances
    )

def test_write_graphviz():
    import boost.graph as bgl
    g = bgl.Graph()
    v1 = g.add_vertex()
    v2 = g.add_vertex()
    e = g.add_edge( v1, v2 )
    v_prop = g.add_vertex_property( 'node_id', 'string' )
    v_prop[ v1 ] = 'hello'
    e_prop = g.add_edge_property( 'node_id', 'string' )
    e_prop[ e ] = 'test'
    g.clear_vertex( v2 )
    g.remove_vertex( v2 )
    g.write_graphviz( 'test.dot' )
test_write_graphviz()
