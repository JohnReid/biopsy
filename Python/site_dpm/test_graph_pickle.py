#
# Copyright John Reid 2008
#

"""
Code to test pickling of boost graphs.
"""

import pickle

def pickle_test(g):
    "Tests pickling by pickling/unpickling a graph"
    print 'Pickling'
    g_pickle_repr = pickle.dumps(g, True)
    print 'Unpickling'
    g_copy = pickle.loads(g_pickle_repr)
    print 'Unpickled'
    assert g_copy.num_vertices() == g.num_vertices()
    assert g.__dict__().keys() == g_copy.__dict__().keys()
    assert g.__dict__().values() == g_copy.__dict__().values()

def test_std_graph():
    import boost.graph as bgl
    g = bgl.Graph()
    #g.add_vertex()
    pickle_test(g)

if '__main__' == __name__:
    test_std_graph()
