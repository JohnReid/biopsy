#
# Copyright John Reid 2006
#

"""
Checks a remome for uniqueness of its sequences.
"""


import biopsy, sys

# filename = sys.argv[1]
filename = 'c:/data/remos/100/100.cleaned'

print 'Loading remome from %s' % ( filename )
r = biopsy.Remome.deserialise( filename )

biopsy.blast_remos( r )
# biopsy.check_remome_uniqueness( r )
