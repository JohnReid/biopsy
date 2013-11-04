#
# Copyright John Reid 2007
#

import _gapped_pssms_hmm as hmm

def weblogo_data_from_dist( dist ):
    """Data for weblogo from pssm distribution"""
    import weblogolib, corebio.seq
    return weblogolib.LogoData.from_counts(
            corebio.seq.unambiguous_dna_alphabet,
            dist * 100
    )

def weblogo_from_dist( dist, filename = 'logo.eps' ):
    """Generate a weblogo from a pssm distribution"""
    import weblogolib
    data = weblogo_data_from_dist( dist )
    options = weblogolib.LogoOptions()
    options.size = weblogolib.LogoSize(
            stack_width = 5.4*12,
            stack_height = 5.4*12*5
    )
    options.color_scheme = weblogolib.std_color_schemes[ "classic" ]
    format = weblogolib.LogoFormat( data, options )
    weblogolib.eps_formatter(
            data,
            format,
            open( filename, 'w' )
    )

def convert_format( source, dest, convert_args = '' ):
    """Converts an image from one format to another

    Uses imagemagick convert program which must be installed

    Args:
            source: input image file
            dest: output image file
    """
    import os
    command = 'convert.exe %s "%s" "%s"' % ( convert_args, source, dest )
    status = os.system( command )
    if status:
        raise RuntimeError(
                'Could not convert "%s" to "%s".\nCommand: %s\nStatus: %d'
                % (
                        source,
                        dest,
                        command,
                        status
                )
        )


def format_weblogo_from_dist( dist, basename, ext, convert_args = '' ):
    """Generate a weblogo from a pssm distribution in format defined by ext"""
    import os.path
    d = os.path.dirname( basename )
    if '' != d and not os.path.exists( d ): os.makedirs( d )
    eps_file = basename + '.eps'
    converted_file = basename + '.' + ext
    if eps_file == converted_file:
        raise RuntimeError( 'Extension should not be same as eps' )
    weblogo_from_dist( dist, eps_file )
    convert_format( eps_file, converted_file, convert_args )
    os.remove( eps_file )


def write_model_svg(
        model,
        name = 'model',
        dir = '.',
        show_rev_comp = False ,
        show_dists = True ,
        edge_lengths = 1.5
):
    print 'Writing model as svg: %s' % name
    p_r_given_pre = model.p_r_given_predecessor

    import boost.graph as bgl
    import os.path, os, numpy

    state_map = hmm.StateMap( model.data.K )
    pre = state_map.predecessor_states()
    g = bgl.Digraph()
    state = g.add_vertex_property( 's', 'integer' )
    label = g.add_vertex_property( 'label', 'string' )
    shapes = g.add_vertex_property( 'shape', 'string' )
    if show_dists: shapefile = g.add_vertex_property( 'shapefile', 'string' )
    style = g.add_vertex_property( 'style', 'string' )
    fillcolor = g.add_vertex_property( 'fillcolor', 'string' )
    vertices = list()
    for s in xrange( state_map.S ):
        # check we want to add this vertex - i.e. check its rev comp attr
        if not show_rev_comp and state_map.c( s ):
            vertices.append( None )
            continue

        # add vertex
        v = g.add_vertex()
        vertices.append( v )
        state[ v ] = s
        label[ v ] = 'k=%d\\ns=%s' % (state_map.k( s ), s)
        if state_map.b( s ): shapes[ v ] = 'rect'
        elif state_map.c( s ): shapes[ v ] = 'diamond'
        else: shapes[ v ] = 'ellipse'
        if state_map.g( s ):
            style[ v ] = 'filled'
            fillcolor[ v ] = 'lightgray'

        # make an image of the distribution for this node
        if show_dists:
            basename = os.path.join( name, 'emission_%d' % s )
            m = state_map.m( s )
            dist = model.var_dist.eta[ m:m+1 ]
            # is it a rev_comp base?
            if state_map.c( s ):
                dist = numpy.array( [ [
                        dist[0][3],
                        dist[0][2],
                        dist[0][1],
                        dist[0][0] ] ] )
            format_weblogo_from_dist( dist, os.path.join( dir, basename ), 'png' )
            shapefile[ v ] = basename + '.png'

    edge_label = g.add_edge_property( 'label', 'string' )
    edge_len = g.add_edge_property( 'len', 'float' )
    for s in xrange( state_map.S ):
        if not show_rev_comp and state_map.c( s ): continue
        for p in pre[s]:
            if not show_rev_comp and state_map.c( int( p ) ): continue
            e = g.add_edge( vertices[p], vertices[s] )
            edge_label[ e ] = '%.3f' % p_r_given_pre[p,s]
            edge_len[ e ] = 1.5
            assert state[ vertices[ p ] ] == p
            assert state[ vertices[ s ] ] == s


    # write as SVG
    dot_filename = '%s.dot' % name
    dot_file = os.path.join( dir, dot_filename )
    svg_filename = '%s.svg' % name
    g.write_graphviz( dot_file )
    label = "\\n".join(
            (
                    "Gaps are gray.",
                    "Reverse complement are diamonds.",
                    "Background state is square",
            )
    )
    wd = os.getcwd()
    os.chdir( dir )
    try:
        os.system(
                "neato %s -Tsvg -o%s -Goverlap=scale \"-Glabel=%s\"" % (
                        dot_filename,
                        svg_filename,
                        label
                )
        )
    finally:
        os.chdir( wd )
