#
# Copyright John Reid 2006
#

import os.path
from env import *

from env import *



_default_db_name_map = {
        'entrezgene': 1,
        'ensembl': 2,
        'unigene': 3,
        'embl': 4,
        'uniprot/sptrembl': 5,
}


_dir = get_aliases_dir()



def _get_dict_from_file( f ):
    """Returns a dictionary mapping the first csv entry to a set of the rest

    I.e. 'key,value1,value2,value3' is converted to

    result[ key ] = set( value1, value2, value3 )

    f can be a string or a filehandle
    """

    # if a string was passed in, convert it to a file handle
    if str == type( f ):
        f = open( f, 'r' )

    result = { }

    for line in f:

        # get the csv fields
        fields = line.strip( '\r\n' ).split( ',' )
        if len( fields ) < 1:
            continue

        # the first is the key
        key = fields[ 0 ]
        result[ key ] = set( fields[ 1: ] )

    return result



def transfac_pssm_2_ensembl_genes():
    """Returns a map from Transfac PSSMs to EnsEMBL genes
    """

    return _get_dict_from_file(
            os.path.join(
                    _dir,
                    'transfac_pssm_2_ensembl_genes.txt' ) )



def transfac_pssm_2_entrez_proteins():
    """Returns a map from Transfac PSSMs to Entrez proteins
    """

    return _get_dict_from_file(
            os.path.join(
                    _dir,
                    'transfac_pssm_2_entrez_proteins.txt' ) )



def ensembl_gene_2_entrez():
    """Returns a map from EnsEMBL genes to Entrez genes
    """

    return _get_dict_from_file(
            os.path.join(
                    _dir,
                    'ensembl_gene_2_entrez.txt' ) )


def ensembl_homologies( key_pred, value_pred = None ):
    """Returns a dictionary mapping genes that match key_pred to those that
    match value_pred.

    If key_pred or value_pred are strings, predicates are formed using
    str.startswith()

    If value_pred is not specified, 'not key_pred' is used
    """

    if str == type( key_pred ):
        key_prefix = key_pred
        key_pred = lambda x: x.startswith( key_prefix )

    if None == value_pred:
        value_pred = lambda x: not key_pred( x )

    if str == type( value_pred ):
        value_prefix = value_pred
        value_pred = lambda x: x.startswith( value_prefix )

    result = { }

    for line in open(
            os.path.join(
                    _dir,
                    'ensembl_homologies.txt' ),
            'r' ):

        fields = line.rstrip( '\r\n' ).split( ',' )
        if len( fields ) < 2:
            continue

        keys = filter( key_pred, fields )
        values = set( filter( value_pred, fields ) )

        for k in keys:
            if result.has_key( k ):
                result[ k ] |= values
            else:
                result[ k ] = values

    return result




def db_name_to_index( db_name ):
    """Implements the default mapping from database names to integer ids
    """

    if str == type( db_name ):
        db_name = db_name.lower()

    if not _default_db_name_map.has_key( db_name ):
        return None
    else:
        return _default_db_name_map[ db_name ]


def get_db_by_name( db_name ):
    db = db_name_to_index( db_name )
    if None == db:
        raise RuntimeError, 'Unknown db name: ' + db_name
    return db



def split_db_ref( db_ref ):
    s = db_ref.split( ':' )
    if len( s ) != 2:
        return ( None, None )
    ( db_name, accession ) = s
    return ( db_name_to_index( db_name ), accession )



class DbAliases:
    """Maps entries in databases to sets of entries in other databases
    """

    def __init__( self ):
        self.aliases = { } # map from from_db to maps


    def get_aliases_for( self, db_from ):

        if str == type( db_from ):
            db = db_name_to_index( db_from )
            if None == db:
                raise RuntimeError, 'Unknown database: ' + db_from
            db_from = db

        if not self.aliases.has_key( db_from ):
            self.aliases[ db_from ] = { }

        return self.aliases[ db_from ]



    def get_map_for( self, db_from, db_to ):

        _map = self.get_aliases_for( db_from )

        if str == type( db_to ):
            db = db_name_to_index( db_to )
            if None == db:
                raise RuntimeError, 'Unknown database: ' + db_to
            db_to = db

        if not _map.has_key( db_to ):
            _map[ db_to ] = { }

        return _map[ db_to ]





    def get( self, db_from, acc_from, db_to ):
        """Get all the accession ids in db_to that are mapped to from acc_from
        """

        _map = self.get_map_for( db_from, db_to )

        if not _map.has_key( acc_from ):
            _map[ acc_from ] = set()

        return _map[ acc_from ]


def get_ensembl_aliases( dbs ):
    """Builds aliases contained in ensembl synonyms file
    """
    _ensembl_synonyms_file = os.path.join(
            biopsy.get_data_dir(),
            'ensembl_synonyms.txt' )

    print 'Building gene ref aliases from', _ensembl_synonyms_file
    result = DbAliases()
    ensembl = get_db_by_name( 'ensembl' )
    if None == ensembl:
        raise RuntimeError, 'Could not find ensembl database id'

    # for each line in the file
    for line in open( _ensembl_synonyms_file, 'r' ):

        # get each csv field
        fields = line.rstrip( '\r\n' ).split( "," )
        if len( fields ):

            # the first is the gene stable id
            gene_stable_id = fields.pop( 0 )
            if len( fields ):

                # the remaining fields are the references
                for ref in fields:

                    ( db, accession ) = split_db_ref( ref )
                    if None == db or not db in dbs:
                        continue
                    # print ref, ( db, accession )

                    # print ( ensembl, gene_stable_id, db )
                    result.get( ensembl, gene_stable_id, db ).add( accession )

    return result


if __name__ == '__main__':
    print transfac_pssm_2_ensembl_genes()[ 'M00025' ]
    print ensembl_gene_2_entrez()[ 'ENSMUSG00000000197' ]
    print ensembl_homologies( 'ENSMUSG' )[ 'ENSMUSG00000000103' ]
    print ensembl_homologies( 'ENSMUSG', 'ENSCAFG' )[ 'ENSMUSG00000000103' ]
