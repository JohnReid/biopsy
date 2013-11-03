#
# Copyright John Reid 2006
#

import biopsy

print biopsy.lookup( 'EntrEz', 'b', 'Entrez' )

def f( arg ):
    print 'In "f"', arg

# t = biopsy.test()
# t.connect_slot( f )
