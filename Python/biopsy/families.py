#
# Copyright John Reid 2006
#

from graph import *

def read_families_file( f ):
    """Reads a paralog file generated from ensembl_homologies.pl

    subtypes can contain a list of paralog subtypes we are interested in
    Yields sequences of genes that form families"""
    if isinstance( f, str ): f = open( f, 'r' )
    for l in f:
        yield l.strip().split(',')

if '__main__' == __name__:
    family_graph = graph_generate(
            read_families_file( 'C:/Data/Ensembl/mouse_families.txt' )
    )
    graph_print_info( family_graph )
    print graph_are_connected(
            family_graph,
            graph_vertex( family_graph, 'ENSMUSG00000000103' ),
            graph_vertex( family_graph, 'ENSMUSG00000049576' )
    )
    print graph_are_connected(
            family_graph,
            graph_vertex( family_graph, 'ENSMUSG00000000001' ),
            graph_vertex( family_graph, 'ENSMUSG00000000149' )
    )
    print graph_are_connected(
            family_graph,
            graph_vertex( family_graph, 'ENSMUSG00000049576' ),
            graph_vertex( family_graph, 'ENSMUSG00000000149' )
    )
