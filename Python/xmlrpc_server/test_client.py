#
# Copyright John Reid 2010
#

"""
Test client for XMLRPC server.
"""

import xmlrpclib, bx.align.maf
from cStringIO import StringIO
from Bio.Seq import Seq, SeqRecord
from Bio.Alphabet import generic_dna


def make_maf(blocks):
    "Add a header to a sequence of blocks to make into MAF format."
    return '##maf version=1\n#\n%s' % '\n\n'.join(blocks)


s = xmlrpclib.ServerProxy('http://localhost:8000')

# Print list of available methods
# print s.system.listMethods()

genome = 'mm9'
chr = 'chr18'
start, end = 79231952, 79232622
blocks = s.maf_extract_range_indexed(genome, '%s.%s' % (genome, chr), start, end)
reader = bx.align.maf.Reader(StringIO(make_maf(blocks)))
alignments = list(reader)
for alignment in alignments:
    #print alignment
    for component in alignment.components:
        segment = '%s:%d,%d' % (component.src, component.start, component.end)
        seq = Seq(component.text, generic_dna)
        print segment
        print seq

    print
