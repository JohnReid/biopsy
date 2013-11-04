#
# Copyright John Reid 2006
#

import numpy, numpy.linalg, math

seqs = [
        'a',
        'c',
        'g',
]

def index_from_nucleo( b ):
    if 'a' == b or 'A' == b: return 0
    if 'c' == b or 'C' == b: return 1
    if 'g' == b or 'G' == b: return 2
    if 't' == b or 'T' == b: return 3
    raise RuntimeError( 'Unknown nucleotide' )

def nucleo_from_index( i ):
    if 0 == i: return 'a'
    if 1 == i: return 'c'
    if 2 == i: return 'g'
    if 3 == i: return 't'
    raise RuntimeError( 'Unknown index' )

def build_counts( seqs, pseudo_counts = 1.0 ):
    l = len(seqs) and len(seqs[0]) or 0
    f = numpy.zeros( (l,4), numpy.float64 )
    f += pseudo_counts
    for s in seqs:
        assert l == len(s)
        for i, b in enumerate( s ):
            f[i,index_from_nucleo(b)] += 1
    return f

pseudo_counts = 1.0
counts = build_counts( seqs, pseudo_counts )
f = counts / ( len( seqs ) + 4.0 * pseudo_counts )
print 'f'
print f

theta = numpy.asarray( [ t for t in f.flat ] )
print 'theta'
print theta

# Fisher information matrix
I = numpy.matrix( numpy.zeros( ( len(theta), len(theta) ), numpy.float64 ) )
for i, t in enumerate( theta ):
    I[i,i] = 1.0 / t
# I inverse
Iinv = numpy.linalg.inv(I)

def V_theta( theta, s ):
    assert len(s) * 4 == len( theta )
    result = numpy.zeros( ( len( theta ), ), numpy.float64 )
    for i, t in enumerate( theta ):
        if nucleo_from_index( i % 4 ) == s[ i / 4 ]:
            result[ i ] = 1.0 / t
    return numpy.asmatrix( result )
print 'V_theta(a)'
print V_theta( theta, 'a' )
print V_theta( theta, 'a' ).shape

def k( s1, s2 ):
    return V_theta( theta, s1 ) * Iinv * V_theta( theta, s2 ).T

print k( 'a', 'c' )
print k( 'a', 'a' )
print k( 'a', 't' )
print k( 't', 't' )
