#
# Copyright John Reid 2006
#

import boost.graph as bgl
import re
import cookbook

go_graph = bgl.Digraph()
id_to_vertex = { }
vertex_to_id = go_graph.add_vertex_property( 'id', 'string' )
def get_vertex( go_id ):
    if not go_id in id_to_vertex:
        v = go_graph.add_vertex()
        id_to_vertex[ go_id ] = v
        vertex_to_id[ v ] = go_id
    return id_to_vertex[ go_id ]

stanza_name_re = re.compile( '^\[(.*)\]$' )
tag_value_re = re.compile( '^([^:]*):\s*([^\{]*)\s*(?:\s*(\S*))' )

def generate_stanzas( f ):
    stanza = cookbook.Bunch( _type_ = 'Header', fields = [] )

    for l in f:
        # remove comments
        i = l.find( '!' )
        if 0 == i: continue
        if -1 != i and l[i-1] != '\\': l = l[:i]

        # remove whitespace at end
        l = l.strip()
        if len( l ) == 0: continue

        # is it a stanza name line?
        match = stanza_name_re.match( l )
        if match:
            yield stanza # yield the last stanza
            stanza = cookbook.Bunch( _type_ = match.group(1), fields = [] )

        # is it a tag value line?
        else:
            match = tag_value_re.match( l )
            if match:
                stanza.fields.append( ( match.group( 1 ), match.group( 2 ) ) )
            else:
                print "Could not parse '%s'" & l

def generate_terms( f ):
    """Generate all the terms in the OBO file"""
    for stanza in generate_stanzas( f ):
        if stanza._type_ != 'Term': continue
        term = cookbook.Bunch(
                id = '',
                is_a = [],
                name = ''
        )
        for f in stanza.fields:
            if 'id' == f[0]: term.id = f[1]
            if 'name' == f[0]: term.name = f[1]
            if 'is_a' == f[0]: term.is_a.append(f[1])
        yield term

for term in generate_terms( open( 'c:\\data\go\\gene_ontology_edit.obo', 'r' ) ):
    v = get_vertex( term.id )
    for v2 in [ get_vertex( is_a ) for is_a in term.is_a ]:
        go_graph.add_edge( v, v2 )

print '# terms', len( go_graph.vertices )
print '# edges', len( go_graph.edges )
