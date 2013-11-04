#
# Copyright John Reid 2006
#
import string
import pickle
import re
import aliases
import os.path
from env import *

_ensembl_synonyms_file = os.path.join( get_biopsy_data_dir(), 'ensembl_synonyms.txt' )
_gene_2_ref_pickle_file = os.path.join( get_biopsy_data_dir(), 'gene_2_ref.map' )
_ref_2_gene_pickle_file = os.path.join( get_biopsy_data_dir(), 'ref_2_gene.map' )
_to_ignore = re.compile( '^AFFY|^AgilentProbe|^IPI|^GO' \
        + '|^RefSeq_peptide_predicted|^RefSeq_dna_predicted' )
_gene_2_ref = None
_ref_2_gene = None

def gene_2_ref_map():
    global _gene_2_ref
    if not _gene_2_ref:
        print 'Loading gene_2_ref_map from', _gene_2_ref_pickle_file
        _gene_2_ref = pickle.load( open( _gene_2_ref_pickle_file, 'r' ) )
    return _gene_2_ref

def ref_2_gene_map():
    global _ref_2_gene
    if not _ref_2_gene:
        print 'Loading ref_2_gene_map from', _ref_2_gene_pickle_file
        _ref_2_gene = pickle.load( open( _ref_2_gene_pickle_file, 'r' ) )
    return _ref_2_gene

def build_gene_ref_aliases():
    """Builds aliases contained in ensembl synonyms file
    """

    print 'Building gene ref aliases from', _ensembl_synonyms_file
    a = aliases.DbAliases()
    ensembl = aliases.get_db_by_name( 'ensembl' )
    if None == ensembl:
        raise RuntimeError, 'Could not find ensembl database'

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

                    ( db, accession ) = aliases.split_db_ref( ref )
                    if None == db:
                        continue
                    # print ref, ( db, accession )

                    # print ( ensembl, gene_stable_id, db )
                    a.get( ensembl, gene_stable_id, db ).add( accession )

    return a


def build_gene_refs():
    """Builds maps from genes to references and back again
    """

    print 'Building gene refs from', _ensembl_synonyms_file
    f = open( _ensembl_synonyms_file, 'r' )
    # print f

    _gene_2_ref = { }
    _ref_2_gene = { }

    # for each line in the file
    for line in f:

        # get each csv field
        fields = line.rstrip( '\r\n' ).split( "," )
        if len( fields ):

            # the first is the gene stable id
            gene_stable_id = fields.pop( 0 )
            if len( fields ):

                # insert a set into the dictionary if necessary
                if not _gene_2_ref.has_key( gene_stable_id ):
                    _gene_2_ref[ gene_stable_id ] = set( fields )
                else:
                    print 'Already seen this gene stable id', gene_stable_id
                    _gene_2_ref[ gene_stable_id ] = _gene_2_ref[ gene_stable_id ] | set( fields )

                # the remaining fields are the references
                for ref in fields:

                    # ignore refs that match our regex
                    if to_ignore.match( ref ):
                        continue

                    # add the set to the dictionary if necessary
                    if not _ref_2_gene.has_key( ref ):
                        _ref_2_gene[ ref ] = set()

                    # add the ref to the set
                    _ref_2_gene[ ref ].add( gene_stable_id )

    # print _gene_2_ref.keys()
    # print _gene_2_ref[ 'ENSMUSG00000000126' ]
    # print _ref_2_gene[ 'Q3UV74' ]

    print 'Found', len( _ref_2_gene ), ' references with genes'
    print 'Found', len( _gene_2_ref ), ' genes with references'

    # find all the refs that have more than one gene id
    if 0:
        print 'Refs with more than one gene associated'
        for k, v in _ref_2_gene.iteritems():
            if len( v ) > 1:
                print k, v

    print 'Pickling maps to', _gene_2_ref_pickle_file, \
            'and', _ref_2_gene_pickle_file
    pickle.dump( _gene_2_ref, open( _gene_2_ref_pickle_file, 'w' ) )
    pickle.dump( _ref_2_gene, open( _ref_2_gene_pickle_file, 'w' ) )
