#
# Copyright John Reid 2010
#

import _biopsy as B

def hit_as_str(hit):
    return '%s %d:%d (%s) %.5f' % (
        hit.binder,
        hit.location.position,
        hit.location.position+hit.location.length,
        hit.location.positive_strand and 'True' or 'False',
        hit.p_binding
    )

hits = B.HitVec()
hits.append(B.Hit('A', B.HitLocation(5,10,True), .5))
hits.append(B.Hit('B', B.HitLocation(5,10,True), .3))
hits.append(B.Hit('C', B.HitLocation(5,10,True), .7))
print map(hit_as_str, hits)

hit_vectors = B.HitsVec()
hit_vectors.append(hits)
hit_vectors.append(hits)

for hits in hit_vectors:
    for hit in hits:
        print hit_as_str(hit)
