#
# Copyright John Reid 2007
#

import biopsy.transfac, csv, cookbook, os

def hit_to_references( hit, predicate ):
    """Takes one hit and yields all the db references associated with it
    that the predicate returns true for"""
    link = biopsy.transfac.TableLink( hit.binder )
    for f in link.entry.factors:
        g = f.link.entry.gene
        if None != g:
            for r in g.entry.db_refs:
                p = predicate( r )
                if p:
                    yield p

def hits_to_identifiers( hits, predicate ):
    """Takes a set of hits and returns a map from ids to hits

    Where the ids are mgi identifiers and the hits are the hits associated with
    them
    """
    result = { }
    for hit in hits:
        for ref in hit_to_references( hit, predicate ):
            if ref in result:
                result[ ref ].add( hit )
            else:
                result[ ref ] = set( [ hit ] )
    return result

def get_ensembl_2_mgi():
    return dict(
            (
                    biopsy.transfac.DbRef.parse_as( row[0], biopsy.transfac.db.ensembl ),
                    [
                            biopsy.transfac.DbRef.parse_as( r, biopsy.transfac.db.mgi )
                            for r in row[1:]
                    ]
            )
            for row in csv.reader(open(os.path.join(biopsy.get_data_dir(), 'identifiers', 'ensembl_mgi_xrefs.csv', "r") ))
    )

def get_entrez_2_mgi():
    # map from mgi ids to entrez
    return cookbook.ReverseDict(
            row[:2]
            for row in csv.reader(open(os.path.join(biopsy.get_data_dir(), 'identifiers', 'entrez_2_mgi.csv')))
    )


class ensembl_mouse_to_mgi( object ):
    """Takes ensembl mouse db refs and converts them to MGI refs"""

    def __init__(self):
        ensembl_2_mgi = get_ensembl_2_mgi()

    def __call__( self, db_ref ):
        if db_ref.db != biopsy.transfac.db.ensembl: return None
        if db_ref.table != 'ENSMUSG': return None
        if db_ref not in self.ensembl_2_mgi: return None
        else: return self.ensembl_2_mgi[ db_ref ][ 0 ].table



def hit_map_to_identifier_map(
        hit_map,
        predicate
):
    """
    Return a dictionary that maps remo names to maps of identifiers to hits
    """
    return dict(
            [
                    (
                            name,
                            hits_to_identifiers(
                                    value,
                                    predicate
                            )
                    )
                    for name, value
                    in hit_map.iteritems()
            ]
    )


def remos_to_identifiers( hit_map ):
    # convert to a map of mgi identifiers
    mgi_identifier_map = hit_map_to_identifier_map(
            hit_map,
            ensembl_mouse_to_mgi()
    )

    # map from mgi ids to entrez
    entrez_2_mgi_map = get_entrez_2_mgi()

    # a map for entrez identifiers
    entrez_identifier_map = dict(
            [
                    (
                            remo,
                            [
                                    entrez_2_mgi_map.reverse[ mgi_id ]
                                    for mgi_id in mgi_ids
                                    if mgi_id in entrez_2_mgi_map.reverse
                            ]
                    )
                    for remo, mgi_ids
                    in mgi_identifier_map.iteritems()
            ]
    )

    return entrez_identifier_map
