#
# Copyright John Reid 2007
#

def graph_as_svg(
        g,
        basename,
        *neato_properties
):
    import os
    dot_filename = '%s.dot' % basename
    svg_filename = '%s.svg' % basename
    g.write_graphviz( dot_filename )
    property_args = " ".join( neato_properties )
    cmd = 'neato %s -Tsvg -o %s %s' % (
            dot_filename,
            svg_filename,
            property_args
    )
    #print 'Running: %s' % cmd
    if os.system( cmd ):
        raise RuntimeError( 'Could not execute: %s' % cmd )
    if os.access(dot_filename, os.R_OK):
        os.remove(dot_filename)
