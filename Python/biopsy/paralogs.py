#
# Copyright John Reid 2006
#

from graph import *

def read_paralog_file( f, subtypes = None ):
    """Reads a paralog file generated from ensembl_homologies.pl

    subtypes can contain a list of paralog subtypes we are interested in
    Yields pairs of genes that are paralogs"""
    if isinstance( f, str ): f = open( f, 'r' )
    counts = { }
    for l in f:
        fields = l.strip().split(',')
        # print fields
        if len(fields) < 3:
            # print 'Short line:', l
            continue
        if not counts.has_key( fields[0] ): counts[ fields[0] ] = 0
        counts[ fields[0] ] += 1
        if ( None != subtypes ) and ( fields[0] not in subtypes ):
            continue # ignore those paralogs we're not interested in
        yield (fields[1], fields[2])
    print
    print "\n".join( [ "%s: %d" % (subtype, c) for subtype, c in counts.iteritems() ] )
    print

if '__main__' == __name__:
    paralog_graph = graph_generate(
            read_paralog_file(
                    'C:/Data/Ensembl/mouse_paralogs.txt',
                    [
                            # 'Fungi/Metazoa group',
                            # 'Coelomata',
                            # 'Bilateria',
                            # 'Chordata',
                            # 'Amniota',
                            # 'Tetrapoda',
                            'Euteleostomi',
                            'Theria',
                            'Eutheria',
                            'Euarchontoglires',
                            'Murinae',
                            'Mus musculus',
                    ]
            )
    )

    graph_print_info( paralog_graph )
    print graph_are_connected(
            paralog_graph,
            graph_vertex( paralog_graph, 'ENSMUSG00000000103' ),
            graph_vertex( paralog_graph, 'ENSMUSG00000049576' )
    )
    print graph_are_connected(
            paralog_graph,
            graph_vertex( paralog_graph, 'ENSMUSG00000000001' ),
            graph_vertex( paralog_graph, 'ENSMUSG00000000149' )
    )
