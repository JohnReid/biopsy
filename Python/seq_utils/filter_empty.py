#!/usr/bin/env python
#
# Copyright John Reid 2010
#

"""
Code that reads in sequences from a FASTA file and outputs the non-empty ones.
"""


import sys, corebio.seq_io.fasta_io as F, corebio.seq, numpy
from optparse import OptionParser


#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
    '-m',
    '--max-sequences',
    dest='max_seqs',
    type='int',
    default=-1,
    help="Set a limit on the number of sequences output."
)
options, args = option_parser.parse_args()


#
# Check args
#
if 1 != len(args):
    print >> sys.stderr, 'USAGE: %s <fasta file>' % __file__
    sys.exit(-1)

fasta = args[0]
if '-' == fasta:
    input = sys.stdin
else:
    input = open(fasta, 'r')

#
# Read the sequences
#
alphabet = corebio.seq.reduced_nucleic_alphabet
for i, seq in enumerate(F.iterseq(input, alphabet)):
    if len(seq):
        F.writeseq(sys.stdout, seq)
