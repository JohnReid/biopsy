#
# Copyright John Reid 2009
#

"""
Some matrices to use as examples for structural basis for gapped PWMs.
"""


import numpy
from itertools import imap, ifilter



cebp_M00912_seqs = """
TTGC-TAAACT
TTGC-GCAAGA
TTGC-TCAATT
TTTC-TCAATA
TTGCTTGAATA
TTGC-CCAACT
"""


mef2_M00941_seqs = """
GGTTATTTTTAA-
GGTTATAATTAA-
AGTTATTTTTAA-
GGTTATTTTTAG-
TGCTAT-TTTAAA
GGCTAT-TTTAAG
ACCTATATTTAG-
AGCTAT-TTTAGG
GGCTATATTTAT-
"""


def parse_site_strings(s):
    """
    Take a string such as:

    GGTTATTTTTAA-
    GGTTATAATTAA-
    AGTTATTTTTAA-
    GGTTATTTTTAG-
    TGCTAT-TTTAAA
    GGCTAT-TTTAAG
    ACCTATATTTAG-
    AGCTAT-TTTAGG
    GGCTATATTTAT-

    and parse it to produce a list of sites.
    """
    return filter(None, imap(str.strip, s.split('\n')))


def char_to_base(c):
    if 'a' == c or 'A' == c:
        return 0
    if 'c' == c or 'C' == c:
        return 1
    if 'g' == c or 'G' == c:
        return 2
    if 't' == c or 'T' == c:
        return 3
    if '-' == c:
        return 4
    raise RuntimeError('Unknown base: %s' % c)


def seqs_to_freqs(seqs):
    if 0 == len(seqs):
        raise RuntimeError('No sequences given')
    K = len(seqs[0])
    counts = numpy.zeros((K, 5))
    for seq in seqs:
        if K != len(seq):
            raise RuntimeError('Not all sequences are the same length')
        for k, b in enumerate(imap(char_to_base, seq)):
            counts[k, b] += 1
    gaps = 1. - counts[:,4] / len(seqs)
    base_counts = counts[:,:4]
    freqs = (base_counts.T / base_counts.sum(axis=1)).T
    return freqs, gaps


def complete_gaps_to_list(gaps):
    "Converts format from array to number for each base to list of tuples (position, gap freq)."
    for i, g in enumerate(gaps):
        if g < 1.:
            yield i, g


if '__main__' == __name__:
    from information_content import MultipleGappedPWM, BG, ic
    for tag, seqs in [
        ('CEBP', cebp_M00912_seqs),
        ('MEF2', mef2_M00941_seqs),
    ]:
        freqs, gaps = seqs_to_freqs(parse_site_strings(seqs))
        pwm = MultipleGappedPWM(freqs, complete_gaps_to_list(gaps))
        logo = pwm.logo()
        logo.save('%s-logo.png' % tag)
        IC = ic(pwm, BG(pwm.K)) / numpy.log(2.)
        print '%s: IC = %f bits' % (tag, IC)
