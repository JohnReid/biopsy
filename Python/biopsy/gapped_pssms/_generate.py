#
# Copyright John Reid 2006
#

import numpy, numpy.random
from _maths import *

def reverse_complement( s ):
    result = numpy.zeros_like( s )
    for i in xrange( len( s ) ):
        result[ len(s) - i - 1 ] = 3 - s[i]
    return result

class GappedPssm( object ):
    def __init__(
            self,
            phi,
            varphi,
            K,
            alpha = [ 1.0, 1.0 ]
    ):
        """Generate a pssm from a base distribution

        K: PSSM length
        phi: dirichlet prior for theta
        varphi: dirichlet prior for pi
        alpha: prior on gamma parameter
        """
        self.theta = numpy.empty( (K+1,4), numpy.float64 )
        for j in xrange( K + 1 ):
            if 0 == j: self.theta[j,:] = sample_from_dirichlet( varphi )
            else: self.theta[j,:] = sample_from_dirichlet( phi )
        self.h = numpy.random.randint(0, K - 1)
        self.gamma = numpy.random.beta( alpha[0], alpha[1] )
        self.K = K

    def __str__( self ):
        return 'Theta:\n%s\nh:%d\ngamma:%f' % (
                str( self.theta ),
                self.h,
                self.gamma
        )

    def sample_from( self ):
        has_gap = bool( numpy.random.random() > self.gamma ) # does it have a gap?
        result = numpy.zeros( self.K + 1, dtype = numpy.int32 )

        # fill in the pssm bases
        for i in xrange( self.K ):
            if has_gap and i > self.h: idx = i + 1
            else: idx = i
            result[ idx ] = sample_from_discrete( self.theta[i+1] )

        # fill in the gap base
        if has_gap: gap_base = self.h + 1
        else: gap_base = self.K
        result[ gap_base ] = sample_from_discrete( self.theta[0] )

        return result, has_gap



def generate_synthetic_sequence( pi, L ):
    """Generate one sequence

    pi: sequence distribution
    L: length
    """
    return numpy.array(
            [
                    sample_from_discrete( pi )
                    for i in xrange( L )
            ],
            dtype = numpy.int32
    )

def base_to_str( base ):
    """Converts 0,1,2,3 to A,C,G,T"""
    if 0 == base: return 'A'
    if 1 == base: return 'C'
    if 2 == base: return 'G'
    if 3 == base: return 'T'
    raise RuntimeError( 'Bad base: %d' % base )

def seq_to_str( seq ):
    """Our sequences are held as arrays, this converts to A,C,G,T strings"""
    return ''.join( [ base_to_str(s) for s in seq ] )

def place_binding_site_in_sequence( seq, pssm ):
    """Replaces part of a sequence with a binding site from the pssm

    seq: sequence
    pssm: pssm

    returns (s,g) where s is the position the site starts at and g is whether
    there is a gap
    """
    sample_seq, has_gap = pssm.sample_from()

    rev_comp = bool( numpy.random.random() > 0.5 ) # is it reverse complemented?
    if rev_comp: sample_seq = reverse_complement( sample_seq )

    s = numpy.random.randint( 0, len( seq ) - pssm.K ) # where in sequence?

    # replace sequence
    seq[s:s+pssm.K+1] = sample_seq

    return s, has_gap, rev_comp

class ModelSample( object ):
    def __init__(
            self,
            phi,
            varphi,
            K,
            N,
            av_length,
            alpha = [ 1.0, 1.0 ],
            verbose = False
    ):
        """Generate some synthetic sequences

        N: number of sequences
        K: length of PSSM
        av_length: expected length
        phi: prior for pssm
        varphi: prior for background
        alpha: prior for gamma
        """
        if verbose:
            print 'Creating PSSM of length %d' % K
        self.pssm = GappedPssm( K = K, phi = phi, varphi = varphi, alpha = alpha )

        if verbose:
            print 'Creating %d sequences of average length %d' % (N, av_length)
        length = 0
        while length < K + 1:
            length = numpy.random.poisson( av_length )
        self.seqs = [
                generate_synthetic_sequence(
                        self.pssm.theta[0,:],
                        length
                )
                for n in xrange( N )
        ]
        self.locations = []
        self.has_gap = []
        self.rev_comp = []
        self.ungapped_sites = []
        for n in xrange( N ):
            s, g, rev_comp = place_binding_site_in_sequence( self.seqs[n], self.pssm )
            self.locations.append( s )
            self.has_gap.append( g )
            self.rev_comp.append( rev_comp )
            # Calculate the pssm that would be generated if the sites were known
            gapped_site = self.seqs[n][s:s+K+1]
            if rev_comp: gapped_site = reverse_complement( gapped_site )
            if g:
                ungapped_site = numpy.concatenate(
                        (
                                gapped_site[:1+self.pssm.h],
                                gapped_site[2+self.pssm.h:]
                        )
                )
            else:
                ungapped_site = gapped_site[:-1]
            # print g, self.pssm.h, gapped_site, ungapped_site
            assert len( ungapped_site ) == self.pssm.K
            self.ungapped_sites.append( ungapped_site )

    def dist_of_sites( self ):
        """The distribution of bases at the actual sites in the sequences

        This will not be the same as the actual pssm
        """
        return dist_from_seqs( self.ungapped_sites )

    def __str__( self ):
        return (
                'Gapped Pssm Model Sample:\n'
                'Pssm: %s\n'
                '# sequences: %d\n'
                'Lengths: %s\n'
                'Starts: %s\n'
                'Has gap: %s'
        ) % (
                str( self.pssm ),
                len( self.seqs ),
                ','.join( [ str(len(s)) for s in self.seqs ] ),
                ','.join( [ str(l) for l in self.locations ] ),
                ','.join( [ str(g) for g in self.has_gap ] ),
        )


def dist_from_seqs( seqs ):
    if not len( seqs ): return numpy.array((), dtype = numpy.float64)
    K = len( seqs[0] )
    result = numpy.zeros( (K,4), dtype = numpy.float64 )
    for i in xrange( K ):
        for s in seqs:
            assert s[i] < 4
            result[i,s[i]] += 1.0
    return result / len( seqs )
