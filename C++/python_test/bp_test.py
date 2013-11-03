#
# Copyright John Reid 2006
#

import sys

print 'Importing bptest'
_path_to_module = '../bin/Bio/Biopsy/_bptest.pyd/vc-8_0/debug/threading-multi'
sys.path.append( _path_to_module )
import _bptest

print 'Creating string vec'
strings = _bptest.StringVec()

print 'Appending string'
strings.append( 'a' )
