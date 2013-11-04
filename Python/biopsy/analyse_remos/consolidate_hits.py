#
# Copyright John Reid 2008
#

"""
Associates BiFa results with their transcription factors and consolidates overlapping
hits.
"""

import biopsy, biopsy.transfac




def map_hits(hits, binder_map):
    """
    Takes a sequence of BiFa binding hits and returns a sequence of BiFa hits over the
    binders that are mapped to binder_map.

    @arg hits: BiFa hits
    @arg binder_map: A dict like object that maps binder names. E.g. mapping BiFa pssm names to transcription factors
    @return: A sequence of BiFa hits
    """
    result = biopsy.HitVec()
    for hit in hits:
        for mapped_binder in binder_map[hit.binder]:
            result.append(biopsy.Hit(mapped_binder, hit.location, hit.p_binding))
    return result




class Pssm2FactorNameMap(dict):
    def __missing__(self, key):
        try:
            link = biopsy.transfac.TableLink(key)
            if biopsy.transfac.trans_data.site == link.db:
                factor_names = [link.entry.factors[0].link.entry.name.upper()] # take first factor name if it is a Site
            else:
                factor_names = [link.entry.factor_name.upper()] # take factor name from Matrix
        except:
            raise RuntimeError("Do not know factor name for binder: %s" % key)
        self[key] = factor_names
        return factor_names




class HitVecDict(dict):
    """
    A dict where missing values are initialised to biopsy.HitVec().
    """
    def __missing__(self, k):
        self[k] = biopsy.HitVec()
        return self[k]




def hits_per_binder(hits):
    """
    Returns a dict mapping binders to HitVecs of the hits for that binder
    """
    result = HitVecDict()
    for hit in hits:
        result[hit.binder].append(hit)
    return result





def hit_location_cmp(hit1, hit2):
    """
    Is the location of hit1 before the location of hit2? Used to sort hits.
    """
    diff = hit1.location.start() - hit2.location.start()
    return 0 != diff and diff or hit1.location.end() - hit2.location.end()





def consolidate_overlapping_hits(hits):
    """
    Yields hits from the sequence of hits where all overlapping hits have been merged into a single hit.
    Ignores the binder value so that the merged hits have the binder of the first overlapping hit.
    The p(binding) of the consolidated hits is the maximum of their constituent hits.
    """
    hits = sorted(hits, cmp=hit_location_cmp) # make sure in location order
    last_hit = None
    for hit in hits:
        if None == last_hit:
            last_hit = hit
        else:
            # does it overlap?
            if last_hit.location.end() > hit.location.start():
                # yes - so update last hit
                last_hit = biopsy.Hit(
                  last_hit.binder,
                  biopsy.HitLocation(
                    last_hit.location.start(),
                    hit.location.end()-last_hit.location.start(),
                    last_hit.location.positive_strand),
                  max(last_hit.p_binding, hit.p_binding)
                )
            else:
                # no - so yield last_hit
                yield last_hit
                last_hit = hit
    if None != last_hit:
        yield last_hit




def consolidate_hits(hits):
    """
    Takes a sequence of BiFa hits and consolidates all of those for the same binder that overlap.

    @arg hits: A sequence of BiFa hits
    @return:  A consolidated sequence of BiFa hits
    """
    from itertools import chain
    result = biopsy.HitVec()
    for hits in hits_per_binder(hits).values():
        result.extend(consolidate_overlapping_hits(hits))
    return result




def maximal_chain_hits(hits, max_boxes=1000000):
    """
    Takes one sequence of hits and returns the maximal chain of non-overlapping hits.
    """
    hit_vec = biopsy.HitsVec()
    hit_vec.append(hits)
    hit_vec.append(hits) # algorithm doesn't work on one sequence of hits, so add it twice
    return biopsy.analyse_max_chain(hit_vec, max_boxes)





if '__main__' == __name__:
    from biopsy import Hit, HitLocation
    def hit_as_str(hit):
        return 'Hit("%s",%.3f,%s)' % (
          hit.binder,
          hit.p_binding,
          hit.location
        )
    Hit.__repr__ = hit_as_str
    def hit_location_as_str(location):
        return 'Loc(%d->%d,%s)' % (
          location.start(),
          location.end(),
          location.positive_strand and '+' or '-'
        )
    HitLocation.__repr__ = hit_location_as_str
    hits = [
      Hit('A', HitLocation(10, 10, True), .5),
      Hit('B', HitLocation(15, 10, True), .6),
      Hit('B', HitLocation(5, 10, True), .4),
      Hit('A', HitLocation(5, 10, True), .7),
    ]
    consolidated = consolidate_hits(hits)
