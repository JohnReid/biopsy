#
# Copyright John Reid 2006-2010
#

import gzip, os.path, re, biopsy

_base_data_dir = os.path.join(biopsy.get_data_dir(), 'ensembl', 'genomes')

def get_genome_dir( genome ):
    return os.path.join( _base_data_dir, genome )

def get_chromosome_file( genome, chromosome ):
    """Returns an open file handle for the chromosome file in the given genome"""
    genome_dir = get_genome_dir( genome )
    files = os.listdir( genome_dir )
    chr_re = re.compile( '\.dna\.chromosome\.%s$' % ( str( chromosome ) ) )
    matched_files = [ f for f in files if chr_re.search( f ) ]
    if 0 == len( matched_files ):
        raise RuntimeError(
                'Did not find chromosome "%s" in directory: %s' % (
                        str( chromosome ),
                        genome_dir
                )
        )
    if len( matched_files ) > 1:
        raise RuntimeError(
                'Expecting only one match in genome directory: %s' % genome_dir
        )
    filename = os.path.join( genome_dir, matched_files[0] )
    print 'File: %s' % filename
    return open( filename, 'r' )

def get_chromosome_sequence( chromosome_file, offset, length ):
    """Gets the sequence starting at the zero-based offset"""
    chromosome_file.seek( offset )
    result = chromosome_file.read( length )
    if len( result ) != length:
        raise RuntimeError(
                'Could not read %d bases from position %d in %s' % (
                        length,
                        offset,
                        str( chromosome_file )
                )
        )
    return result

if '__main__' == __name__:
    chromosome_file = get_chromosome_file( 'mus_musculus_core_28_33d', 6 )
    print get_chromosome_sequence( chromosome_file, 6589024, 77 )
