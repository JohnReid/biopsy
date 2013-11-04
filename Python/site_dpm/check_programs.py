#
# Copyright John Reid 2009
#

"""
Code to execute to check all is ok with programs.
"""


from itertools import chain, imap


def hits_for_gene(ucsc_analysis, gene):
    "Yield all the hits for a gene in the analysis."
    return chain(*(ucsc_analysis[gene]))


def binder_for_hit(hit):
    "@return: The hit's binder."
    return hit.binder


def mapper(map):
    "@return: A function that maps keys through the given map."
    def map_fn(key):
        "@return: The value for the key."
        return map[key]
    return map_fn


if '__main__' == __name__:
    import ucsc_promoters, analysis, biopsy, shared

    tag = 'lower'
    try: ucsc_analysis
    except NameError: ucsc_analysis = ucsc_promoters.get_analysis(tag)
    pssm_map = analysis.get_pssm_to_ensembl_map_min_range()
    ensembl_names = shared.get_all_ensembl_names()
    for gene in [
        'ENSMUSG00000000581',
        'ENSMUSG00000035403',
        'ENSMUSG00000022184',
        'ENSMUSG00000061099',
        'ENSMUSG00000031949',
        'ENSMUSG00000037447',
        'ENSMUSG00000026753',
        'ENSMUSG00000022877',
        'ENSMUSG00000066273',
        'ENSMUSG00000031131',
        'ENSMUSG00000015305',
        'ENSMUSG00000061911',
    ]:
        print '%s : %s' % (gene, ','.join(map(mapper(ensembl_names), imap(mapper(pssm_map), imap(binder_for_hit, hits_for_gene(ucsc_analysis, gene))))))
