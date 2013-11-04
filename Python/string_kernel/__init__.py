
from _string_kernel import *

def convert_to_sequence( s ):
    seq = Sequence()
    for x in s:
        seq.append( x )
    return seq

def convert_to_sequence_vec( s ):
    seqs = SequenceVec()
    for x in s:
        seqs.append( convert_to_sequence( x ) )
    return seqs

def normalise( k ):
    import numpy
    """Normalises a string kernel"""
    k_copy = numpy.array( k, copy = True )
    sqrt = numpy.sqrt( k.diagonal() )
    for i in xrange( k.shape[0] ):
        for j in xrange( k.shape[1] ):
            k_copy[i,j] = k[i,j] / ( sqrt[i] * sqrt[j] )
    return k_copy
