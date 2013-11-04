#
# Copyright John Reid 2006
#

import elementtree.ElementTree as ET
import graph
import entrez
import os.path

"""Code to handle Psi-Mi protein-protein interaction files
"""

interaction_path = \
        './entry/{net:sf:psidev:mi}interactionList/{net:sf:psidev:mi}interaction'

interactor_path = \
        '{net:sf:psidev:mi}participantList/*/*/{net:sf:psidev:mi}xref/' \
        + '{net:sf:psidev:mi}primaryRef'

negative_path = \
        './entry/{net:sf:psidev:mi}interactionList/{net:sf:psidev:mi}interaction/' \
        + '{net:sf:psidev:mi}negative'

def parse( handle ):
    """Parses an input stream to produce a list of interactions.

    Each interaction is itself a list of database references
    """

    # parse the input
    tree = ET.parse( handle )

    # for each interaction
    num_interactions = 0
    for interaction in tree.getroot().findall( interaction_path ):
        interaction_list = []
        for interactor in interaction.findall( interactor_path ):
            interaction_list.append(
                    interactor.attrib.get( "db" )
                    + ':'
                    + interactor.attrib.get( "id" ) )
        yield interaction_list
        num_interactions += 1
    print '# interactions: %d' % num_interactions

    # for each negative value
    parent_map = dict((c, p) for p in tree.getiterator() for c in p)
    for neg in tree.getroot().findall( negative_path ):
        if 'true' == neg.text:
            print 'Negative:', parent_map[ neg ].attrib[ 'id' ]



def convert_to_edges( interactions ):
    """Takes a list of interactions and converts them to edges

    Where edges = pairwise interaction
    """
    # for each interaction
    num_edges = 0
    for interaction in interactions:

        # loop through all pairs of interactors
        for ( i, interactor1 ) in enumerate( interaction ):
            for interactor2 in interaction[ i + 1: ]:
                yield ( interactor1, interactor2 )
                num_edges += 1
    print '%d edges' % num_edges


def duplicate_vertices( edges ):
    """Returns a new edge list with 2 copies of each node

    Inserts 4 edges, one for each connection between the nodes

    We may need two copies of each node to handle the case where a protein
    interacts with itself
    """

    for e in edges:
        v0a = 'A_' + e[ 0 ]
        v0b = 'B_' + e[ 0 ]
        v1a = 'A_' + e[ 1 ]
        v1b = 'B_' + e[ 1 ]
        # print v0a, v0b, v1a, v1b
        yield ( v0a, v1a )
        yield ( v0b, v1a )
        yield ( v0a, v1b )
        yield ( v0b, v1b )



def split_ref( ref ):
    """Splits a database reference into a db name and an accession number
    """
    ( db_name, acc_id ) = ref.split( ':' )
    return ( db_name, acc_id )



def strip_edges( edges, db_name ):
    """Only retains edges that link refs in the given database and don't have
    -1 as an accession number
    """

    dbs = set()
    num_edges = 0
    for ( p1, p2 ) in edges:
        db1, acc1 = split_ref( p1 )
        db2, acc2 = split_ref( p2 )
        if db_name == db1 and db_name == db2 and "-1" != acc1 and "-1" != acc2:
            yield ( acc1, acc2 )
            num_edges += 1
    print '%d edges between proteins in %s' % (num_edges, db_name)



def build_graph( f ):
    print 'Building graph'
    g = graph.build(
            strip_edges(
                    convert_to_edges(
                            parse( f )
                    ),
                    "Entrez Protein"
            )
    )

    # gene_map = g.vertex_properties[ 'genes' ] = g.vertex_property_map( 'string' )
    # entrez.build_vertex_2_gene_map( g, g.vertex_properties[ 'node_id' ], gene_map )

    return g
