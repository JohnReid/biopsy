#
# Copyright John Reid 2006
#

from _biopsy import *

def parse_fasta( f ):
    """Parses a fasta file

    f can either be a file name or a file handle.

    Returns a dictionary mapping keys to sequences

    Not optimised for large files
    """

    if str == type( f ):
        f = open( f, 'r' )

    results = { }
    key = None
    for line in f:
        line = line.rstrip( '\r\n' )
        if line.startswith( '>' ):
            key = line.lstrip( '> ' )
            results[ key ] = ''
        else:
            results[ key ] += line

    return results


def get_char_for_nucleo( nucleotide ):
    return nucleotide
