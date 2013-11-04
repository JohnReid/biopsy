#
# Copyright John Reid 2006
#

import elementtree.ElementTree as ET
import Bio.EUtils as EUtils
import Bio.EUtils.HistoryClient as HistoryClient

def qualifiers( node ):
    "Yields all the qualifiers beneath the node as name, value tuples"
    for qual in node.findall( './/GBQualifier' ):
        name = qual.find( 'GBQualifier_name' )
        if None == name: continue
        value = qual.find( 'GBQualifier_value' )
        if None == value: continue
        yield ( name.text, value.text )

_client = None
def client():
    "Returns a singleton HistoryClient we can work with"
    global _client
    if None == _client:
        _client = HistoryClient.HistoryClient()
    return _client

def filter_db_xrefs( iterable, xref_db = "MGI" ):
    "Yields only those db_xrefs that are for the given db"
    for n, v in iterable:
        if 'db_xref' == n:
            fields = v.split( ':' )
            if 2 != len( fields ): continue
            if xref_db == fields[0]:
                yield fields[1]

def seq_ids( node ):
    "Gets all the seq ids in the given node"
    for seq_id in node.findall( './/GBSeqid' ):
        yield seq_id.text.split( '|' )

def gi_id( node ):
    "Gets the GI id for the node"
    for seq_id in seq_ids( node ):
        if 'gi' == seq_id[0]: return seq_id[1]

def post( ids, entrez_db = "protein" ):
    "Make a request to NCBI, yield each response separately"
    result = client().post( EUtils.DBIds("protein", ids) )
    tree = ET.parse( result.efetch() )
    for seq_node in tree.findall( './/GBSeq' ):
        yield seq_node

def xrefs( ids, entrez_db = "protein", xref_db = "MGI" ):
    "Yield tuples (id, (xrefs))"
    for result in post( ids, entrez_db = entrez_db ):
        yield (
                gi_id( result ),
                [
                        xref
                        for xref
                        in filter_db_xrefs(
                                qualifiers( result ),
                                xref_db = xref_db
                        )
                ]
        )
