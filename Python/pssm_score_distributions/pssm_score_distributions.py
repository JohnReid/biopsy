#
# Copyright John Reid 2007
#

import biopsy, itertools

num_buckets = 100
num_pssms = 5
num_sequences = 200

def pssm_accs():
    "The pssms we will test"
    pssm_filter = biopsy.transfac.PssmFilter(
            #name_regex_pattern = 'NFAT_Q4_01'
    )
    return biopsy.get_transfac_pssm_accessions( pssm_filter )

def sequences_from_remome( remome ):
    "The sequences we will use to test the score distributions"
    for aligned in remome.get_aligned_sequences():
        for remo in remome.get_remos_for( aligned ):
            for species in remo.get_sequence_ids():
                yield remo.get_sequence_for( species, True )

def histogram( acc, score_counts ):
    import pylab, numpy
    pylab.clf()
    pylab.bar(
            xrange( num_buckets ),
            numpy.power( score_counts, 0.25 )
    )
    pylab.savefig( 'graphs/%s.png' % acc )

try:
    remome
except NameError:
    remome_file = 'c:/data/remos/100/100.filtered'
    print 'Loading remome: %s' % remome_file
    remome = biopsy.Remome.load( remome_file )

score_counts = { }
for acc in itertools.islice( pssm_accs(), num_pssms ):
    score_counts[acc] = numpy.zeros( num_buckets, numpy.uint8 )
    for seq in itertools.islice( sequences_from_remome( remome ), num_sequences ):
        hits = biopsy.HitVec()
        p_bind = biopsy.score_pssm_on_sequence(
                pssm_name = acc,
                threshold = 0.0,
                sequence = seq,
                result = hits
        )
        for h in hits:
            score_counts[acc][ int( h.p_binding * num_buckets ) ] += 1
    histogram( acc, score_counts[acc] )
