#
# Copyright John Reid 2013
#

"""
Code to deal with MEME and MEME file formats.
"""

import biopsy
from Bio import SeqIO
from itertools import chain
from collections import defaultdict


def name_matcher(name):
    "Create a function that matches strings in lower case"
    name = name.lower()
    def matcher(other):
        return -1 != other.lower().find(name)
    return matcher


def match_factor(matcher, factor):
    "Does the factor's name or synonyms match the matcher function"
    if matcher(factor.name):
        return True
    else:
        for synonym in factor.synonyms:
            if matcher(synonym):
                return True


def find_matrices(name):
    matcher = name_matcher(name)
    for matrix in biopsy.transfac.Matrix.all():
        for facc in matrix.factors:
            factor = facc.link.entry
            if factor.gene is not None and 'MOUSE' == factor.gene.entry.species:
                if match_factor(matcher, factor):
                    yield matrix, factor
                    break



def logo(dist, tag, dir):
    "Generate a logo with the given tag in the given directory."
    import weblogolib as W
    import corebio.seq as S
    data = W.LogoData.from_counts(S.unambiguous_dna_alphabet, dist)
    options = W.LogoOptions(
        logo_title=tag,
        color_scheme=W.colorscheme.nucleotide,
        show_xaxis=False,
        show_yaxis=True,
        show_fineprint=False,
    )
    format = W.LogoFormat(data, options)
    filename = 'logo-%s' % tag
    #W.eps_formatter(data, format, open(os.path.join(dir, '%s.eps' % filename), 'w'))
    W.png_formatter(data, format, open(os.path.join(dir, '%s.png' % filename), 'w'))
    #W.pdf_formatter(data, format, open(os.path.join(dir, '%s.pdf' % filename), 'w'))


def dist_for_pssm(pssm):
    "@return: The PSSM's frequencies."
    import numpy as N
    return N.array(
        [
            [pssm.dists[i].get_freq(b) for b in xrange(4)]
            for i in xrange(len(pssm.dists))
        ]
    )



def look_for_matrices(names):
    for name in names:
        print name
        for matrix, factor in find_matrices(name):
            print matrix.acc, matrix.name, factor.acc, factor.name
            logo(dist_for_pssm(biopsy.get_pssm(str(matrix.acc))), '%s-%s' % (name, matrix.acc), 'logos')


def write_minimal_meme_matrix(out, acc):
    """
    The minimal MEME format for a motif looks something like::

        MOTIF crp
        letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009
        0.000000  0.176471  0.000000  0.823529
        0.000000  0.058824  0.647059  0.294118
        0.000000  0.058824  0.000000  0.941176
        0.176471  0.000000  0.764706  0.058824
        0.823529  0.058824  0.000000  0.117647
        0.294118  0.176471  0.176471  0.352941
        0.294118  0.352941  0.235294  0.117647
        0.117647  0.235294  0.352941  0.294118
        0.529412  0.000000  0.176471  0.294118
        0.058824  0.235294  0.588235  0.117647
        0.176471  0.235294  0.294118  0.294118
        0.000000  0.058824  0.117647  0.823529
        0.058824  0.882353  0.000000  0.058824
        0.764706  0.000000  0.176471  0.058824
        0.058824  0.882353  0.000000  0.058824
        0.823529  0.058824  0.058824  0.058824
        0.176471  0.411765  0.058824  0.352941
        0.411765  0.000000  0.000000  0.588235
        0.352941  0.058824  0.000000  0.588235
    """
    pssm_info = biopsy.get_pssm(acc)
    print >>out, (
        "MOTIF %s %s\n"
        "letter-probability matrix: alength= 4 w= %d nsites= %d E= %e\n"
        "%s\n"
    ) % (
        biopsy.get_pssm_name(acc), acc,
        len(pssm_info.dists), pssm_info.sites, 0.,
        "\n".join(
            '  '.join(("%.6f" % dist.get_freq(b)) for b in xrange(4))
            for dist in pssm_info.dists
        )
    )


def write_minimal_meme(out, accs, bg_freqs=None):
    """
    The header looks something like::

        MEME version 4

        ALPHABET= ACGT

        strands: + -

        Background letter frequencies
        A 0.303 C 0.183 G 0.209 T 0.306
    """
    if bg_freqs is None:
        bg_freqs = [.25, .25, .25, .25]

    print >>out, (
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: + -\n\n"
        "Background letter frequencies\n"
        "A %.3f C %.3f G %.3f T %.3f \n\n"
    ) % (
        bg_freqs[0], bg_freqs[1], bg_freqs[2], bg_freqs[3],
    )
    for acc in accs:
        write_minimal_meme_matrix(out, acc)
        print >>out, ""


