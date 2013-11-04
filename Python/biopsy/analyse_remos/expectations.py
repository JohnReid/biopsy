#
# Copyright John Reid 2008
#

"""
Calculates expectations of the number of hits on a per binder basis
"""

class _FloatDict(dict):
    def __missing__(self, k):
        self[k] = 0.0
        return self[k]

class _IntDict(dict):
    def __missing__(self, k):
        self[k] = int(0)
        return self[k]

def expected_num_hits_per_binder(hits, expectations=None):
    """
    Calculates expectations of the number of hits on a per binder basis.

    @return: dict mapping binder names to expectations
    """
    if None == expectations:
        expectations = _FloatDict()
    for hit in hits:
        expectations[hit.binder] += hit.p_binding
    return expectations

def num_hits_per_binder(hits, counts=None):
    """
    Calculates expectations of the number of hits on a per binder basis.

    @return: dict mapping binder names to counts
    """
    if None == counts:
        counts = _IntDict()
    for hit in hits:
        counts[hit.binder] += 1
    return counts
