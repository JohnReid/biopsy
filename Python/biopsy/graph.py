#
# Copyright John Reid 2006
#

import tempfile
import os.path

def build( edges ):
    import boost.graph as bgl

    # Create a graph to build
    g = bgl.Graph()

    # Construct some property maps
    node_ids = g.add_vertex_property( 'node_id', 'string' )
    weights = g.add_edge_property( 'weight', 'float' )

    # create the edges and vertices
    name_2_vertex = { }
    for ( name1, name2 ) in edges:

        # have we seen name 1 before?
        if not name_2_vertex.has_key( name1 ):
            v1 = g.add_vertex()
            name_2_vertex[ name1 ] = v1
            node_ids[ v1 ] = name1
        else:
            v1 = name_2_vertex[ name1 ]

        # have we seen name 2 before?
        if not name_2_vertex.has_key( name2 ):
            v2 = g.add_vertex()
            name_2_vertex[ name2 ] = v2
            node_ids[ v2 ] = name2
        else:
            v2 = name_2_vertex[ name2 ]

        # is the edge in the graph already?
        e = g.edge( v1, v2 )
        if None == e:
            e = g.add_edge( v1, v2 )
            weights[ e ] = 1
        else:
            weights[ e ] += 1

    return g


def create_vertex_index_map( g ):
    idx = 0
    indices = g.vertex_properties[ 'vertex_indices' ] = g.add_vertex_property( 'integer' )
    for v in g.vertices:
        indices[ v ] = idx
        idx = idx + 1
    return indices

def print_vector_in_hb_format(f, vec):
    num_on_line = 0
    for entry in vec:
        f.write( '%5d' % ( entry, ) )
        num_on_line = num_on_line + 1
        if num_on_line == 16:
            f.write('\n')
            num_on_line = 0
    if num_on_line != 0:
        # pad to end
        f.write(' '*((16 - num_on_line)*5))
        f.write('\n')



def write_harwell_boeing_format(
        g,
        f,
        title = '<No title>',
        key = '<No key>' ):

    # our format only handles numbers up to 5 digits
    if g.num_vertices() > 99999:
        raise NameError, 'Too many vertices'

    # names = g.get_vertex_object_map("vertex_names")

    # for each vertex - calculate the column indices and row indices
    colidx = 1
    colptr = []
    rowind = []
    colptr.append(colidx)
    vertex_indices = create_vertex_index_map( g )
    for v in g.vertices:
        # print names[v]
        for a in g.adjacent_vertices(v):
            if vertex_indices[a] <= vertex_indices[v]:
                # print ' :', names[a]
                colidx = colidx + 1
                rowind.append(vertex_indices[a] + 1)
        colptr.append(colidx)

    # line 1
    f.write('%72.72s%8.8s\n' % (title, key))

    # line 2
    ptrcrd = (len(colptr) * 5) / 80 + 1
    indcrd = (len(rowind) * 5) / 80 + 1
    totcrd = ptrcrd + indcrd
    valcrd = rhscrd = 0
    f.write(
            '%14d%14d%14d%14d%14d\n'
            % (totcrd, ptrcrd, indcrd, valcrd, rhscrd))

    # line 3
    mxtype = 'PSA'
    nrow = ncol = g.num_vertices()
    nnzero = g.num_edges()
    neltvl = 0
    f.write(
            '%3.3s%11.11s%14d%14d%14d%14d\n'
            % (mxtype, '', nrow, ncol, nnzero, neltvl))

    # line 4
    ptrfmt = indfmt = '(16I5)'
    valfmt = rhsfmt = ''
    f.write(
            '%16.16s%16.16s%16.16s%16.16s\n'
            % (ptrfmt, indfmt, valfmt, rhsfmt))

    # the column pointers
    print_vector_in_hb_format( f, colptr )

    # the row indices
    print_vector_in_hb_format( f, rowind )

    return vertex_indices




def write_matrix_market_format( g, f, predicate = lambda v: True ):
    """Writes the graph, g, in matrix market format to filehandle, f.

    Only writes those vertices for which predicate is true.
    Returns ( num_vertices, num_non_zeros, vertex_index_dict )
    """

    # do we have weights associated with the edges?
    weights = None
    if g.edge_properties.has_key( 'weight' ):
        weights = g.edge_properties[ 'weight' ]

    if None != weights:
        f.write( '%%MatrixMarket matrix coordinate real symmetric\n' )
    else:
        f.write( '%%MatrixMarket matrix coordinate pattern symmetric\n' )

    # get a dictionary from vertices to indices     and calc num_vertices
    vertex_indices = g.add_vertex_property( 'integer' )
    index_vertices = { }
    num_vertices = 0
    for v in filter( predicate, g.vertices ):
        num_vertices += 1
        vertex_indices[ v ] = num_vertices
        index_vertices[ num_vertices ] = v

    # calculate num non zeros in matrix
    edge_set = set()
    for v in filter( predicate, g.vertices ):
        for a in filter( predicate, g.adjacent_vertices( v ) ):
            if vertex_indices[ a ] <= vertex_indices[ v ]:
                edge_set.add( ( vertex_indices[ a ], vertex_indices[ v ] ) )
    num_non_zeros = len( edge_set )

    f.write(
            str( num_vertices )
            + ' '
            + str( num_vertices )
            + ' '
            + str( num_non_zeros )
            + '\n' )

    for e in edge_set:
        line = \
                str( e[ 0 ] ) \
                + ' ' \
                + str( e[ 1 ] )
        if None != weights:
            line += \
                    ' ' \
                    + str( weights[ \
                            g.edge( \
                                    index_vertices[ e[ 0 ] ], \
                                    index_vertices[ e[ 1 ] ] ) ] )
        f.write( line + '\n' )

    return ( num_vertices, num_non_zeros, vertex_indices, index_vertices )



def write_as_connected_components( g, file_prefix ):
    import boost.graph as bgl
    components = g.add_vertex_property( 'integer' )
    num_components = bgl.connected_components( g, components )
    for i in range( num_components ):
        filename = file_prefix + '_' + str( i ) + '.mm'
        f = open( filename, 'w' )
        # print f
        ( num_vertices, num_non_zeros, vertex_indices ) = \
                write_matrix_market_format(
                        g,
                        f,
                        lambda v: components[ v ] == i )
        # print i, num_vertices, num_non_zeros
    return ( components, num_components )


def get_sub_graph_copy( g_orig, vertex_pred, edge_pred = lambda x: True ):
    """Copy those vertices and edges that satisfy the predicates
    """
    import boost.graph as bgl

    f = tempfile.TemporaryFile()
    g_orig.write_graphviz( f )
    f.rewind()
    g = bgl.Graph.read_graphviz( f )
    f.close()

    for e in filter( edge_pred, g.edges ):
        g.remove_edge( e )

    for v in filter( edge_pred, g.vertices ):
        g.clear( v )
        g.remove( v )

    return g


def graph_2_distance_map( g, components_dir, name, beta, eps ):
    """Calculates the distances between all nodes in the graph

    Returns a dict where tuples of node_ids in the graph are the keys and
    the distances are the values. Distances below eps are not returned.

    beta is the parameter for the diffusion kernel on the graph

    components_dir is a directory to place working files in

    name is a string to use when generating file names
    """
    import boost.graph as bgl

    #
    # Calculate the connected components
    #
    # print 'Calculating connected components'
    components = g.add_vertex_property( 'integer' )
    num_components = bgl.connected_components( g, components )
    # print num_components, 'connected component(s)'


    result = { }

    #
    # For each component
    #
    if not os.path.exists( components_dir ):
        os.mkdir( components_dir )
    components_prefix = os.path.join( components_dir, name )
    for i in range( num_components ): # for each component

        # for testing just do one component
        # if 2 != i: continue


        # create an adjacency matrix for this component
        adj_filename = components_prefix + '_' + str( i ) + '_adj.mm'
        predicate = lambda v: components[ v ] == i
        ( num_vertices, num_non_zeros, vertex_indices, index_vertices ) = \
                write_matrix_market_format(
                        g,
                        open( adj_filename, 'w' ),
                        predicate )


        # we need a map from adjacency matrix row/col indices to nodes
        adjacency_nodes = { }
        idx = 1
        for v in filter( predicate, g.vertices ):
            adjacency_nodes[ idx ] = v
            idx += 1

        # create and run matlab script to calculate the kernel
        kernel_filename = components_prefix + '_' + str( i ) + '_ker.mm'
        if not os.path.exists( kernel_filename ):
            tmp_script = "\
A = mmread( '%ADJACENCY_FILE%' );\n\
size( A )\n\
L = calc_laplacian( A );\n\
K = calc_kernel( L, %BETA% );\n\
K_sp = sparsify( K, %EPSILON% );\n\
mmwrite( '%KERNEL_FILE%', K_sp );\n\
exit\n"
            script_text = tmp_script. \
                    replace( '%ADJACENCY_FILE%', adj_filename ). \
                    replace( '%EPSILON%', str( eps ) ). \
                    replace( '%BETA%', str( beta ) ). \
                    replace( '%KERNEL_FILE%', kernel_filename )
            script_filename = 'tmp_script.m'
            open( script_filename, 'w' ).write( script_text )
            matlab_command = 'C:\\apps\\MATLAB7\\bin\\win32\\MATLAB.exe /nosplash ' \
                    + '/r ' + script_filename.strip( '.m' )
            # print matlab_command
            os.system( matlab_command )



        # once we have the kernel we can populate the edge weights
        # are the kernel's values
        f = open( kernel_filename, 'r' )
        n = m = nnz = 0
        node_ids = g.vertex_properties[ 'node_id' ]
        for line in f: # this for loop reads a matrix market format file

            if line.startswith( '%' ): continue

            line = line.strip( '\n' )
            split = line.split()
            # Is this the first line of data with matrix sizes?
            if n == 0:
                ( n, m, nnz ) = [ int( x ) for x in split ]
                # print n, m, nnz
            else:
                idx1 =  int( split[ 0 ] )
                idx2 =  int( split[ 1 ] )
                prot1 = node_ids[ adjacency_nodes[ idx1 ] ]
                prot2 = node_ids[ adjacency_nodes[ idx2 ] ]

                # we are only interested in 'A' -> 'B' weights
                if prot1.startswith( 'A_' ) and prot2.startswith( 'B_' ):
                    prot1 = prot1[ 2: ]
                    prot2 = prot2[ 2: ]
                    distance = float( split[ 2 ] )

                    # only put value in if prot1 < prot2
                    if prot1 <= prot2:
                        result[ ( prot1, prot2 ) ] = distance

    return result

def graph_has_vertex( g, i ):
    """Returns True if graph has vertex with this id"""
    return i in g.id_to_vertex

def graph_vertex( g, i, add_if_necessary = False ):
    """Returns vertex with given id, i"""
    if add_if_necessary and i not in g.id_to_vertex:
        v = g.add_vertex()
        g.id_to_vertex[ i ] = v
        g.vertex_properties[ 'vertex_id' ][ v ] = i
    return g.id_to_vertex[ i ]

def graph_print_info( g ):
    print '# vertices: %d' % len( g.vertices )
    print '# edges: %d' % len( g.edges )
    print '# connected components: %d' % graph_component_map( g )[1]

def graph_component_map( g ):
    """Returns (component_map, num_components)
    """
    import boost.graph as bgl
    try: g.component_map
    except:
        g.component_map = g.add_vertex_property( 'components', 'integer' )
        g.num_components = bgl.connected_components( g, g.component_map )
    return ( g.component_map, g.num_components )

def graph_are_connected( g, v1, v2 ):
    """Are vertices, v1 and v2, connected in graph g?

    Does lazy evaluation of connected components and stores result in graph
    """
    component_map, num_components = graph_component_map( g )
    return component_map[ v1 ] == component_map[ v2 ]

def graph_generate( ids ):
    """Creates a graph from a sequence of sequences

    Expects a sequence of sequences of ids.

    Returns a graph with ids as nodes and connections between successive ids
    in each sequence

    The graph will have an id_to_vertex dictionary attribute
    """
    import boost.graph as bgl
    g = bgl.Graph()
    g.id_to_vertex = { }
    g.add_vertex_property( 'vertex_id', 'string' )
    for seq in ids: # join the sequence up in a line
        last_v = None
        for i in seq:
            v = graph_vertex( g, i, add_if_necessary = True )
            if None != last_v and last_v not in g.adjacent_vertices( v ):
                g.add_edge( v, last_v )
            last_v = v
    return g
