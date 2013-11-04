#
# Copyright John Reid 2008
#

"""
Code to read TRANSFAC sites from file.
"""

import os
from cookbook.named_tuple import namedtuple

def read_sites(f):
    seqs = []
    matrix = None
    name = None
    for line in f:
        line = line.strip()
        if not line:
            if seqs:
                yield matrix, name, seqs
            seqs = []
            matrix = None
            name = None
        else:
            if not matrix:
                matrix, name = line.split('\t')
            else:
                seqs.append(line)

def convert_base(b):
    if 'a' == b or 'A' == b: return 0
    if 'c' == b or 'C' == b: return 1
    if 'g' == b or 'G' == b: return 2
    if 't' == b or 'T' == b: return 3
    raise RuntimeError('Unknown base: %s' % str(b))

def convert_seq(seq):
    return map(convert_base, seq)

sites_filename = os.path.join('c:\\', 'Dev', 'MyProjects', 'Vincent', 'transfac_matrix_sites.txt')

TransfacSiteSet = namedtuple('TransfacSiteSet', 'matrix name seqs')
def load_transfac_sites(sites_filename=sites_filename):
    return dict(
      (matrix, TransfacSiteSet(matrix=matrix, name=name, seqs=seqs))
      for matrix, name, seqs
      in read_sites(open(sites_filename))
    )


def convert_transfac_sites(sites):
    result = dict()
    for matrix, site_set in sites.iteritems():
        try:
            converted_seqs = map(convert_seq, site_set.seqs) # convert to feature values
            if not converted_seqs: # if no sequences ignore
                continue
            for s in converted_seqs: # check all same length
                if len(s) != len(converted_seqs[0]):
                    print '%s: not all sequences same length' % matrix
                    break
            else:
                result[matrix] = converted_seqs
        except:
            print 'Problem converting sequences for', matrix
    return result


if '__main__' == __name__:
    import pylab

    sites = list(read_sites(open(sites_filename)))
    for matrix, name, seqs in sites:
        print matrix
        print name
        print sites
        print
        break
    pylab.hist([len(seqs) for matrix, name, seqs in sites], bins=20)
