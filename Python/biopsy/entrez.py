#
# Copyright John Reid 2006
#

from urllib2 import urlopen, URLError
from xml.parsers.expat import ExpatError
import xml.etree.cElementTree as ET
import sys, re, biopsy, generefs
from cookbook import lru_cache

_url_base = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/'

def _get_and_parse_url( url ):
    tree = 0
    while not tree:
        try:
            tree = ET.parse( urlopen( url ) )
        except URLError:
            # print "Retrying...", url
            pass
    return tree

def proteins_from_gene( entrez_gene_acc ):
    result = []
    url = _url_base + 'elink.fcgi?dbFrom=gene&db=protein&cmd=gene_protein&id=' \
            + str( entrez_gene_acc )
    tree = _get_and_parse_url( url )
    for element in tree.findall( './LinkSet/LinkSetDb/Link/Id' ):
        result.append( int( element.text ) )
    return result


def genes_from_protein( entrez_protein_acc ):
    result = []
    url = _url_base + 'elink.fcgi?dbFrom=protein&db=gene&cmd=protein_gene&id=' \
            + str( entrez_protein_acc )
    tree = _get_and_parse_url( url )
    for element in tree.findall( './LinkSet/LinkSetDb/Link/Id' ):
        result.append( int( element.text ) )
    return result


def uniprot_from_gene( entrez_gene_acc ):
    result = []
    # print 'gene entrez acc =', entrez_gene_acc
    url = _url_base + 'efetch.fcgi?db=gene&retmode=xml&id=' \
            + str( entrez_gene_acc )
    # sys.stdout.write( "url = \"" + url + "\"\n" )

    tree = _get_and_parse_url( url )

    parent_map = dict( (c, p) for p in tree.getiterator() for c in p )
    for element in tree.getiterator( 'Dbtag_db' ):
        if 'UniProt' == element.text:
            parent = parent_map[ element ]
            for value in parent.getiterator( 'Object-id_str' ):
                # print value.text
                result.append( value.text )
    return result

def uniprot_from_protein( entrez_protein_acc ):
    result = []
    for gene in genes_from_protein( entrez_protein_acc ):
        result.extend( uniprot_from_gene( gene ) )
    return result

def build_vertex_2_gene_map( g, node_ids, gene_map ):
    ref_2_gene = generefs.ref_2_gene_map()
    num_found = 0
    num_not_found = 0
    for v in g.vertices:
        genes = []
        acc = node_ids[ v ]
        found = 0
        for entrez_gene in genes_from_protein( acc ):
            str_ref = 'EntrezGene:' + str( entrez_gene )
            if ref_2_gene.has_key( str_ref ):
                for gene in ref_2_gene[ str_ref ]:
                    genes.append( gene )
                    found = 1
        if found:
            num_found = num_found + 1
        else:
            num_not_found = num_not_found + 1
        gene_map[ v ] = ','.join( genes )
    print '# found', num_found
    print '# not found', num_not_found


def get_protein_url( entrez_protein_acc ):
    return 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=' \
            + str( entrez_protein_acc )

_name_re = re.compile( '(.*) \\[.*' )
@lru_cache(cache_storage_file=biopsy.serialisation_dir()+'entrez_protein_name.cache')
def protein_name( entrez_protein_acc ):
    """Returns the gene name from the protein accession"""
    global _protein_name_cache
    import xml.etree.ElementTree as ET
    import Bio.EUtils.ThinClient
    eutils = Bio.EUtils.ThinClient.ThinClient()
    dbids = Bio.EUtils.DBIds( "protein", [ str( entrez_protein_acc ) ] )
    tree = ET.parse( eutils.efetch_using_dbids( dbids, retmode='xml', rettype='gp' ) )
    for qual in tree.findall( './/GBQualifier' ):
        name = qual.find( 'GBQualifier_name' )
        if None != name and name.text == 'gene':
            value = qual.find( 'GBQualifier_value' )
            if None != value:
                return value.text
    return None

if __name__ == '__main__':
    print proteins_from_gene( 72058 )

    # for gene in get_genes_from_entrez_protein( 4557669 ):
    #       print get_uniprot_refs_from_entrez_gene( gene )
