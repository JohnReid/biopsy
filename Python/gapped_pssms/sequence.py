#
# Copyright John Reid 2008
#

"""
Code to load and work with sequences.

Sequences can be in two forms, string or numpy:

string: 'acggt'
numpy: numpy.array([0,1,2,2,3], dtype=int)
"""

import numpy
from sequence import numpy

def base_to_index(base):
    "@return: The given base ('a', 'c', 'g', 't', or 'n') as an index [0,4]."
    if 'a' == base or 'A' == base: return 0
    if 'c' == base or 'C' == base: return 1
    if 'g' == base or 'G' == base: return 2
    if 't' == base or 'T' == base: return 3
    if 'n' == base or 'N' == base: return 4
    raise RuntimeError('Unknown base: "%s"' % base)

def index_to_base(base):
    "@return: The given base given as an index [0,4] to a character ('a', 'c', 'g', 't', or 'n')."
    if 0 == base: return 'a'
    if 1 == base: return 'c'
    if 2 == base: return 'g'
    if 3 == base: return 't'
    if 4 == base: return 'n'
    raise RuntimeError('Unknown base: %d' % base)

def seq_to_numpy(seq):
    "Converts a sequence held as a string to numpy array form."
    result = numpy.empty(len(seq), dtype=numpy.uint32)
    for i, s in enumerate(seq):
        result[i] = base_to_index(s)
    return result

def numpy_to_seq(seq):
    "Converts a sequence held as a numpy array to string form."
    return ''.join(index_to_base(b) for b in seq)

def complement(b):
    "@return: The complement of the given base."
    return (4 == b) and 4 or (3 - b)

def rev_comp(seq):
    "@return: The reverse complement of the given sequence (in numpy form)."
    result = numpy.empty_like(seq)
    for i, b in enumerate(seq):
        result[-i-1] = complement(b)
    return result

def num_bases(seqs):
    "@return: The number of bases in the sequences."
    return sum(len(seq) for seq in seqs)

def num_known_bases(seqs):
    "@return: The number of known bases in the sequences."
    return sum(sum( seq != 4 ) for seq in seqs)

def random_sequence(length):
    "@return: a random sequence of the given length."
    import numpy.random
    return numpy.random.randint(4, size=length)

def random_sequences(num_seqs, average_length):
    """
    @arg num_seqs: The number of sequences to return.
    @arg average_length: The expected length of each sequence. (the lengths are actually distributed according to a Poisson distribution).
    @return: A list of random sequences.
    """
    import numpy.random
    return [
      random_sequence(length=numpy.random.poisson(lam = average_length))
      for i in xrange(num_seqs)
    ]

def strip_ns(seq):
    "@return: The sequence stripped of leading and trailing 'N's."
    return seq.strip('nN')

def yield_fasta_sequences(filename):
    "@return: Yield the sequences from a FASTA file."
    from corebio.seq_io.fasta_io import iterseq, dna_alphabet
    for s in iterseq(open(filename, 'r'), dna_alphabet):
        yield s

def has_length(s):
    "@return: True iff length != 0"
    return len(s) != 0

def load_fasta_sequences(filename):
    "@return: A list of sequences, each stripped of leading and trailing 'N's."
    return filter(has_length, map(strip_ns, yield_fasta_sequences(filename)))

def convert_fasta_sequences(filename):
    "@return: A list of sequences, each stripped of leading and trailing 'N's and converted to numpy form."
    return map(seq_to_numpy, load_fasta_sequences(filename))
