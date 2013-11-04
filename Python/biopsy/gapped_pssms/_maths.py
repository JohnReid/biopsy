#
# Copyright John Reid 2006
#

import numpy

def binary_entropy( p ):
    """Entropy of bernoulli trial"""
    return - p * math.log( p ) - ( 1 - p ) * math.log( 1 - p )

def discrete_entropy( P ):
    """Returns entropy of discrete distribution given by P"""
    return sum (
            [
                    - p * math.log( p )
                    for p in P
                    if p != 0.0
            ]
    )

def dirichlet_entropy( U ):
    """Returns entropy of dirichlet distribution"""
    u0 = sum(U)
    return - sum(
            [
                    math.log( scipy.special.gamma( u ) )
                    for u in U
            ]
    ) + math.log(
            scipy.special.gamma( u0 )
    ) + sum(
            [
                    (u - 1.0) * scipy.special.digamma( u )
                    for u in U
            ]
    ) - ( u0 - len(U) ) * scipy.special.digamma( u0 )

def finite_and_a_number( arg ):
    """Returns true iff arg is finite and is not an nan

    Works with numpy arrays
    """
    return numpy.isfinite( arg ).all() and not numpy.isnan( arg ).any()

def sample_from_dirichlet( alpha ):
    """Sample from dirichlet distribution

    See - http://en.wikipedia.org/wiki/Dirichlet_distribution

    alpha: parameters
    """
    y = numpy.array( [
            numpy.random.gamma( a )
            for a in alpha
    ] )
    s = sum( y )
    y /= s
    return y

def expected_entropy_under_dirichlet( U ):
    """http://www.carleton-scientific.com/isipta/PDF/021.pdf - eq(5)"""
    s = numpy.sum( U )
    psi = scipy.special.digamma( 1.0 + s )
    return sum(
            [
                    u / s * ( psi - scipy.special.digamma(u + 1.0) )
                    for u in U
            ]
    )

def test_expected_entropy_under_dirichlet( u, num_samples = 1000 ):
    h = 0.0
    for i in xrange( num_samples ):
        h += discrete_entropy( sample_from_dirichlet( u ) )
    h /= num_samples
    print h
    print expected_entropy_under_dirichlet( u )
#test_expected_entropy_under_dirichlet( .01 * numpy.ones( 4, numpy.float64 ) )

def sample_from_discrete( P ):
    """Samples from a discrete distribution

    Returns 0,...,n-1 where n = len(p)

    P: probabilities of each of 0..n-1 outcomes
    """
    x = numpy.random.random()
    for i, p in enumerate( P ):
        if p > x: return i
        x -= p
    raise RuntimeError( 'P sums to less than 1' )

def probabilities_from_logs( logs ):
    """Takes logarithms of probabilities (up to a constant) and returns the
    probabilities.
    """
    # make the largest 0 - as we can get underflow problems
    largest = max( logs )
    e = numpy.exp( logs - largest )
    total = numpy.sum( e )
    return e / total
