#
# Copyright John Reid 2006
#

from _biopsy import *

def _find_sequence_in_remome(
        remome,
        sequence,
        masked = True,
        examine_all_species = False
):
    """Finds a given sequence in the remome

    Returns a list of (aligned_sequence_set, remo) tuples.
    """
    if str == type( remome ):
        remome = Remome.deserialise( remome )
    for aligned in remome.get_aligned_sequences():
        # print aligned.centre_sequence
        for remo in remome.get_remos_for( aligned ):
            for species in remo.get_sequence_ids():
                if examine_all_species or aligned.centre_sequence == species:
                    s = remo.get_sequence_for( species, masked )
                    if -1 != s.find( sequence ):
                        yield ( aligned, remo )
                        break
Remome.find_sequence = _find_sequence_in_remome
