
import os, sys, cookbook
from env import *

_ensembl_data_dir = os.path.join(get_data_dir(), 'ensembl')

def gene_location_filename( genome ):
    return os.path.normpath(
            os.path.join(
                    _ensembl_data_dir,
                    'gene_locations-%s.txt' % genome
            )
    )

def ensembl_src_dir( version ):
    return 'C:\\Dev\\ThirdParty\\perl\\ensembl\\ensembl-%d' % version

def set_perl_lib_path( version ):
    perl5lib = ';'.join(
            [
                    'C:\\Dev\\ThirdParty\\perl\\ensembl\\src\\bioperl-live',
                    '%s\\ensembl\\modules' % ensembl_src_dir( version ),
                    '%s\\ensembl-compara\\modules' % ensembl_src_dir( version ),
                    '%s\\ensembl-variation\\modules' % ensembl_src_dir( version ),
            ]
    )
    os.putenv( 'PERL5LIB', perl5lib )

def call_gene_locations_perl( genome ):
    "Must be in directory where gene_locations.pl is to call this"
    filename = 'gene_locations-%s.txt' % genome
    command = 'perl gene_locations.pl %s >%s' % ( genome, filename )
    print "Command: '%s'" % command
    return os.system( command )

default_versions = [
        ( 28, 'mus_musculus_core_28_33d' ),
        ( 31, 'mus_musculus_core_31_33g' ),
        ( 32, 'mus_musculus_core_32_34' ),
]

def get_gene_locations_from_ensembl( versions ):
    for version, genome in versions:
        set_perl_lib_path( version )
        call_gene_locations_perl( genome )

def gene_locations( genome ):
    "Returns a dict mapping gene ids to bunch(chr, start, end, strand)"
    result = { }
    for l in open( gene_location_filename( genome ), 'r' ):
        fields = l.strip().split(':')
        if len( fields ) != 5:
            print 'Bad record: %s' % l.strip()
        else:
            result[fields[0]] = cookbook.Bunch(
                    chromosome = fields[1],
                    start = int( fields[2] ),
                    end = int( fields[3] ),
                    positive_sense = int( fields[4] ) == 1
            )
    return result
