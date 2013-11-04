#
# Copyright John Reid 2006
#

"""Loads paths between different transfac pssms from a file and returns them in
a dict
"""

import os


class PssmPaths:
    _pssm_path_filename = 'pssm_paths.txt'
    _pssm_paths = None

    @classmethod
    def get_pssm_path_file( cl ):
        return os.path.join( get_aliases_dir(), cl._pssm_path_filename )

    @classmethod
    def read_pssm_paths( cl, f ):
        """Returns a map from tuples of pssms to paths"""
        if str == type( f ):
            f = open( f, 'r' )
        result = { }
        for line in f:
            fields = line.rstrip( '\r\n' ).split( "," )
            if len( fields ) == 0:
                continue
            if len( fields ) < 4:
                raise RuntimeError, 'Too short a line: ' + line
            pssm_1 = fields.pop( 0 )
            pssm_2 = fields.pop( 0 )
            result[ ( pssm_1, pssm_2 ) ] = fields
        return result

    @classmethod
    def get( cl ):
        """Returns a map from pssm tuples to shortest paths"""
        if None == cl._pssm_paths:
            cl._pssm_paths = cl.read_pssm_paths( cl.get_pssm_path_file() )
        return cl._pssm_paths
