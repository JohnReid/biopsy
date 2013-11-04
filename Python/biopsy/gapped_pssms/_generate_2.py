#
# Copyright John Reid 2006
#

import numpy, numpy.random
from _maths import *
from _generate import *


class Site( object ):
    """Represents a site's location, orientation and length in a sequence"""
    def __init__( self, start, rev_comp, length ):
        self.start = start
        self.rev_comp = rev_comp
        self.length = length

    def to_seq_coord( self, model_idx ):
        if self.rev_comp: return self.start + self.length - model_idx - 1
        else: return self.start + model_idx

    def to_model_idx( self, seq_coord ):
        if self.rev_comp: return self.start + self.length - seq_coord - 1
        else: return seq_coord - self.start


class GappedSite( Site ):
    "Represents a site of a gapped pssm"
    def __init__( self, start, rev_comp, length, gap_idx ):
        Site.__init__( self, start, rev_comp, length )
        self.gap_idx = gap_idx

    def r( self, idx ):
        """Returns:
                -1 if idx is outside the pssm
                0 if idx hits the gap
                [1, self.length) otherwise"""
        if idx < 0 or idx >= self.length: return -1 # are we outside the pssm?
        if idx == self.gap_idx + 1: return 0 # did we hit the gap?
        if idx <= self.gap_idx: return idx + 1 # are we below the gap?
        else: return idx # we are above the gap

    def idx( self, r ):
        """Returns index into the site for given r"""
        assert( 0 <= r )
        assert( self.length > r )
        if 0 == r: return self.gap_idx + 1 # did we hit the gap?
        if r - 1 <= self.gap_idx: return r - 1 # are we below the gap?
        else: return r # we are above the gap




def random_site( seq_length, model_length ):
    return Site(
            start = numpy.randint( seq_length - model_length ),
            rev_comp = numpy.uniform() > 0.5,
            length = model_length
    )


def generate_sequences( num_seqs, pssm, av_length ):
    "Yields sequences with sites embedded"
    for i in xrange( num_seqs ):
        # first generate a sequence without a site
        seq = generate_synthetic_sequence(
                pssm.theta[0,:],
                numpy.random.poisson( av_length )
        )
        # pick somewhere to put the site
        site = random_site( len( seq_without_site ) )
        # change the bases in the sequence where the site is
        for i in xrange( pssm.K ):
            idx = s + pssm.offset( i, g, rev_comp )
            sample = sample_from_discrete( pssm.theta[i+1] )
            if rev_comp: seq[ idx ] = 3 - sample
            else: seq[ idx ] = sample
