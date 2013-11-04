#
# Copyright John Reid 2007
#

import biopsy.pssm_paths, biopsy.entrez

try:
    paths = biopsy.pssm_paths.PssmPaths.get()
except:
    import warnings
    warnings.warn('Could not get shortest protein paths between PSSMs')

def shortest_paths( seq_1, seq_2, max_length ):
    "Yields the shortest protein interaction paths between the 2 pssms."
    already_examined = set()
    for pssm_1 in seq_1:
        for pssm_2 in seq_2:
            if pssm_1 < pssm_2:
                pair = ( pssm_1, pssm_2 )
            else:
                pair = ( pssm_2, pssm_1 )
            if pair in already_examined: continue
            else: already_examined.add( pair )
            if paths.has_key( pair ):
                path = paths[ pair ]
                if len( path ) <= max_length:
                    yield pair, path


def path_graph(
        binders_1,
        binders_2,
        max_length
):
    """Returns a boost.graph representing the protein links between the lists of hits"""
    import boost.graph as bgl
    g = bgl.Graph()
    id_map = { }
    g.add_vertex_property( name = 'id', type = 'string' )
    g.add_vertex_property( name = 'label', type = 'string' )
    g.add_vertex_property( name = 'fillcolor', type = 'color' )
    g.add_vertex_property( name = 'style', type = 'string' )
    g.add_vertex_property( name = 'shape', type = 'string' )
    g.add_vertex_property( name = 'URL', type = 'string' )
    already_examined = set()
    for b1 in binders_1:
        for b2 in binders_2:

            # the pair needs to have smaller binder first
            reversed = b1 >= b2
            if not reversed: pair = ( b1, b2 )
            else: pair = ( b2, b1 )

            # only look at pairs we have not seen already
            if pair in already_examined: continue
            else: already_examined.add( pair )

            # look for path between binders
            if not paths.has_key( pair ): continue
            path = paths[ pair ]

            # check max length
            if len( path ) > max_length: continue

            # iterate over path adding vertices and edges
            last_protein = None
            for protein in path:
                if not( protein in id_map ):
                    v = g.add_vertex()
                    id_map[protein] = v
                    g.vertex_properties['id'][v] = protein

                    url = biopsy.entrez.get_protein_url( protein )
                    if url: g.vertex_properties['URL'][v] = url

                    label = biopsy.entrez.protein_name( protein )
                    g.vertex_properties['label'][v] = label or 'Unknown'
                else:
                    v = id_map[protein]
                if None != last_protein:
                    last_vertex = id_map[last_protein]
                    if None == g.edge( v, last_vertex ):
                        g.add_edge( v, last_vertex )
                last_protein = protein

            # set colours for first and last in path
            first_v = id_map[ path[ 0 ] ]
            last_v = id_map[ path[ -1 ] ]
            if reversed: first_v, last_v = last_v, first_v
            g.vertex_properties['fillcolor'][ first_v ] = bgl.Color.gray
            g.vertex_properties['style'][ first_v ] = 'filled'
            g.vertex_properties['shape'][ last_v ] = 'diamond'
            print biopsy.transfac.TableLink( b2 ).name, g.vertex_properties['label'][first_v]
            print biopsy.transfac.TableLink( b1 ).name, g.vertex_properties['label'][last_v]
            print


    return g

class PathGraphBuilder( object ):
    """Builds a boost.graph representing paths between hits."""

    def __init__( self ):
        import boost.graph as bgl
        import os
        import biopsy
        self.protein_vertex_map = { }
        self.g = bgl.Graph()
        self.g.add_vertex_property( name = 'id', type = 'string' )
        self.g.add_vertex_property( name = 'label', type = 'string' )
        self.g.add_vertex_property( name = 'fillcolor', type = 'color' )
        self.g.add_vertex_property( name = 'style', type = 'string' )
        self.g.add_vertex_property( name = 'shape', type = 'string' )
        self.g.add_vertex_property( name = 'URL', type = 'string' )
        self.entrez_map = biopsy.aliases._get_dict_from_file(
                os.path.join( biopsy.get_aliases_dir(), 'transfac_pssm_2_entrez_proteins.txt' )
        )

    def add_hits( self, binders, name ):
        for b in binders:
            v = self.get_vertex( b )

    def add_protein( self, protein ):
        """Adds protein and paths to graph. Returns vertex for protein"""
        v = get_vertex( protein )
        if v: return v # already added


    def get_vertex( self, protein ):
        if id in self.protein_vertex_map: return self.protein_vertex_map[ id ]
        else: return None
