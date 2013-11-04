#
# Copyright John Reid 2006
#

import biopsy

for factor in biopsy.get_transfac_factors():
    fragments = biopsy.get_fragments_for_factor( factor )
    if len(fragments):
        print '%s has %d fragments' % (
                factor,
                len( fragments )
        )
        fasta = open( 'fragments/%s.fa' % factor, 'w' )
        for frag in fragments:
            fasta.write( '> %s\n' % frag )
            fasta.write( '%s\n' % biopsy.get_fragment_sequence( frag ) )
        fasta.close()
