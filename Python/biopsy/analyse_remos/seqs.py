#
# Copyright John Reid 2007
#

import biopsy

def convert_dict_seqs_for_phylo_analysis(
        d,
        centre_key = '',
        centre_key_regex = None
):
    """
    Takes a dictionary mapping keys to sequences and returns sequence of sequences, centre first.
    """
    centre_keys = filter(centre_key_regex.match, d.keys())
    if len(centre_keys) > 1:
        raise RuntimeError('Too many centre sequences found: %s' % "; ".join(centre_keys))
    if len(centre_keys) == 0:
        raise RuntimeError('No sequence matches the centre key regex')
    centre_key = centre_keys[0]
    sequences = biopsy.SequenceVec()
    sequences.append(d[centre_key]) # add the centre sequence first
    sequences.extend(s for (k, s) in d.iteritems() if k != centre_key) # add all the other sequences
    assert len(d) == len(sequences), 'Should have same number of sequences in input and output'
    return sequences
