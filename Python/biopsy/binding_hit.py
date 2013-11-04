#
# Copyright John Reid 2006
#

from _biopsy import *


def _hit_str( hit ):
    return ",".join( [
            hit.binder,
            str( hit.location.position ),
            str( hit.location.positive_strand ),
            str( hit.p_binding )
    ] )
Hit.__str__ = _hit_str


def _location_start( location ):
    return location.position
HitLocation.start = _location_start

def _location_end( location ):
    return location.position + location.length
HitLocation.end = _location_end

def _location_overlap( location1, location2 ):
    """Do two hits overlap?"""
    if location1.position < location2.position:
        return location1.end() > location2.position
    else:
        return location2.end() > location1.position
HitLocation.overlap = _location_overlap

def _location_separation( location1, location2 ):
    """The separation between two locations"""
    if location1.position >= location2.end():
        return location1.position - location2.end()
    else:
        return location2.position - location1.end()
HitLocation.separation = _location_separation


def _hits_str( hits ):
    return '\n'.join( [ str( hit ) for hit in hits ] )
HitVec.__str__ = _hits_str



def get_char_for_hit( hit ):
    return hit.binder

def get_score_for_hit( hit ):
    # return math.log( hit.p_binding )
    return hit.p_binding


def get_max_p_binding_over_hits( hits ):
    """Takes a list of hits and returns a dictionary mapping binder names to max( p(binding) ) across all hits"""
    result = { }
    for hit in hits:
        if not result.has_key( hit.binder ) or result[hit.binder] < hit.p_binding:
            result[hit.binder] = hit.p_binding
    return result

def find_pair_in_analysis(
        analysis,
        pair,
        max_separation = None,
        separation = None
):
    """Finds in which analyses a pair of TFs bind

    analysis: Analysis
    pair: A tuple ( binder1, binder2, orientation1, orientation2 )
    max_separation: If specified determines maximum separation
    separation: If specified determines exact separation (over-rides max_separation)

    Returns a list of keys for the analyses
    """
    result = { }
    for k in analysis.get_keys():
        hits = analysis.get_hits_for( k )
        found_pairs = find_pair_in_hits( hits, pair, max_separation, separation )
        if found_pairs:
            result[ k ] = found_pairs
    return result




def find_pair_in_hits(
        hits,
        pair,
        max_separation = None,
        separation = None
):
    """Finds the locations where a pair of TFs bind in a sequence of hits

    hits: The hits
    pair: A tuple ( binder1, binder2, orientation1, orientation2 )
    max_separation: If specified determines maximum separation
    separation: If specified determines exact separation (overrides max_separation)

    returns a sequence of pairs of hits that satisfy the criteria
    """
    ( binder1, binder2, orientation1, orientation2 ) = pair
    result = [ ]
    for h1 in hits:
        if binder1 != h1.binder: continue
        for h2 in hits:
            if binder2 != h2.binder: continue
            if h1.location.overlap( h2.location ): continue
            distance = h1.location.separation( h2.location )
            if None != separation and separation != distance: continue
            if None != max_separation and max_separation < distance: continue
            if h1.location.position < h2.location.position:
                if (
                        h1.location.positive_strand != orientation1
                        or
                        h2.location.positive_strand != orientation2
                ): continue
            else:
                if (
                        h1.location.positive_strand == orientation1
                        or
                        h2.location.positive_strand == orientation2
                ): continue
            result.append( ( h1, h2 ) )
    return result


def hit_over_threshold_predicate(threshold):
    "@return: A function that returns True if the hit is over the threshold given."
    def predicate(hit):
        "@return: True iff the hit's score is above the threshold."
        return hit.p_binding >= threshold
    return predicate

def hits_above_threshold(hits, threshold):
    "@return: Those hits above the threshold."
    return filter(hit_over_threshold_predicate(threshold), hits)
