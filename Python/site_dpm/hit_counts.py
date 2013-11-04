#
# Copyright John Reid 2009
#

"""
Code to convert binding site analyses into hit counts ready for the HDPM.
"""

from shared import *
import biopsy, analysis as A
from itertools import imap
import biopsy.analyse_remos.expectations as EX


class HitAboveThreshold(object):
    "Callable that tests if hit's score is above a given threshold."
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, hit):
        return hit.p_binding >= self.threshold


def hit_mapper(map):
    def map_hit(hit):
        return biopsy.Hit(
            map[hit.binder],
            biopsy.HitLocation(
                hit.location.position,
                hit.location.length,
                hit.location.positive_strand
            ),
            hit.p_binding
        )
    return map_hit

def map_binders(hits, map):
    result = biopsy.HitVec()
    result.extend(imap(hit_mapper(map), hits))
    return result

def count_hits(hit_lists):
    hit_counts = None
    for hits in hit_lists:
        if hits:
            hit_counts = EX.num_hits_per_binder(hits, hit_counts)
    return hit_counts

def convert_analysis_to_counts(analysis, pssm_map, threshold=0.0):
    for key, hit_lists in analysis.iteritems():
        if len(hit_lists):
            yield (
                key,
                count_hits(
                    map_binders(biopsy.hits_above_threshold(hits, threshold), pssm_map)
                    for hits in hit_lists
                )
            )

def remove_empty_counts(hit_counts):
    to_del = set()
    for key, counts in hit_counts.iteritems():
        if None == counts:
            to_del.add(key)
    for key in to_del:
        del hit_counts[key]


@log_exceptions()
@caching_decorator('hit-counts')
def get_hit_counts():
    if options.use_ucsc_seqs:
        import ucsc_promoters
        analysis = ucsc_promoters.get_analysis()
    else:
        raise RuntimeError('Do not know which sequences to use.')
    hit_counts = dict()
    logging.info('Using binding site threshold of %f to calculate hit counts.', options.hit_count_threshold)
    hit_counts.update(convert_analysis_to_counts(analysis, A.get_pssm_to_ensembl_map_min_range(), options.hit_count_threshold))
    remove_empty_counts(hit_counts)
    return hit_counts.copy()


if '__main__' == __name__:
    get_hit_counts()
