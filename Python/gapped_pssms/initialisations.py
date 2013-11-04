#
# Copyright John Reid 2008
#

"""
Examines which of a set of different K-mers to initialise a gapped PSSM with might be best.
"""

import logging, time
from gapped_pssms.sequence import numpy_to_seq, seq_to_numpy, rev_comp


def hamming_distance(k_mer_1, k_mer_2):
    "@return: The Hamming distance between the 2 K-mers (of the same length)."
    if len(k_mer_1) != len(k_mer_2):
        raise ValueError("Can only measure Hamming distance between K-mers of same length")
    return sum(k_mer_1 != k_mer_2)



class K_mer_distance(object):
    """
    A distance metric between 2 K-mers. Allows for shifting and reverse complementing.
    """

    def __init__(self, allowed_shifts=2, shift_cost=1):
        """
        @arg allowed_shifts: How many times we can shift one K-mer while comparing it to the other."
        @arg shift_cost: The edit distance of a shift."
        """

        self.allowed_shifts = allowed_shifts
        "How many times we can shift one K-mer while comparing it to the other."

        self.shift_cost = shift_cost
        "The edit distance of a shift."

    def _distance_wo_shifting(self, k_mer_1, k_mer_2):
        "@return: The distance without considering shifts."
        return min(
          hamming_distance(k_mer_1, k_mer_2),
          hamming_distance(k_mer_1, rev_comp(k_mer_2))
        )

    def __call__(self, k_mer_1, k_mer_2):
        "@return: The distance between the K-mers."
        return min(
          self._distance_wo_shifting(k_mer_1, k_mer_2), # unshifted
          min(
            i * self.shift_cost + self._distance_wo_shifting(k_mer_1[i:], k_mer_2[:-i]) # shift one way
            for i in xrange(1, self.allowed_shifts+1)
          ),
          min(
            i * self.shift_cost + self._distance_wo_shifting(k_mer_1[:-i], k_mer_2[i:]) # shift other way
            for i in xrange(1, self.allowed_shifts+1)
          ),
        )


class DistanceFilter(object):
    """
    Filters K-mers based on how close they are to K-mers that have already passed this filter.
    """

    def __init__(self, distance, min_distance=2):
        """
        @arg distance: A distance measure.
        @arg min_distance: Minimum distance k-mer must be away from K-mers that have already passed.
        """
        self.distance = distance
        "A distance measure."

        self.min_distance = min_distance
        "Minimum distance k-mer must be away from already evaluated K-mers."

        self.passed = set()
        "Those K-mers that have already passed the filter."


    def _distance_to_passed(self, k_mer):
        return min(self.distance(k_mer, seq_to_numpy(other_k_mer)) for other_k_mer in self.passed)

    def __call__(self, k_mer):
        """
        @return: True if and only if the K-mer passes the filter.
        """
        if not self.passed:
            result = True
        else:
            d = self._distance_to_passed(k_mer)
            result = d >= self.min_distance
            if result:
                logging.debug('%s is %d away from previous k-mers', numpy_to_seq(k_mer), d)
        if result == True:
            self.passed.add(numpy_to_seq(k_mer))
        return result




def yield_k_mers(sequences, K):
    """
    @return: Yield the (K-mer, count) in order such that the mers with highest number of occurences come first.
    """
    from hmm import ReverseComplementCollapsingCounter, count_mers
    from heapq import heapify, heappop
    import time

    # find all K-mers collapsed with their reverse complements
    logging.info('Finding all %d-mers in sequences', K)
    start = time.time()
    nmer_counts = ReverseComplementCollapsingCounter(K)
    count_mers(sequences, n=K, callback=nmer_counts)
    logging.info('Took %f seconds to find %d-mers', time.time()-start, K)


    start = time.time()
    counts = list((-count, i, K_mer) for i, (K_mer, count) in enumerate(nmer_counts.counts()))
    heapify(counts)
    logging.info('Took %f seconds to heapify', time.time()-start)
    #import IPython; IPython.Debugger.Pdb().set_trace()

    while counts:
        count, i, K_mer = heappop(counts)
        yield K_mer, -count


def most_frequent_k_mers(sequences, K, max_to_return):
    result = list()
    for K_mer, count in yield_k_mers(sequences, K):
        if 0 == max_to_return:
            break
        result.append((K_mer, count))
        max_to_return -= 1
    return result


if '__main__' == __name__:
    #
    # Test distance metric between K-mers
    #
    from gapped_pssms.sequence import seq_to_numpy
    distance = K_mer_distance()
    K_mers = [
      ('acgtacgtg', 'gacgtacgt', 1), # 1 shift
      ('acgtacgtg', 'ggacgtacg', 2), # 2 shifts
      ('ggacgtacg', 'acgtacgtg', 2), # 2 shifts other way
      ('aaaaaaaaa', 'aaaaaaaaa', 0), # identical
      ('aaaaaaaaa', 'ttttttttt', 0), # reverse complement
      ('acgtacgtc', 'aaaaaaaaa', 7), # very different
    ]
    for k1, k2, d in K_mers:
        n1 = seq_to_numpy(k1)
        n2 = seq_to_numpy(k2)
        assert d == distance(n1, n2) # make sure is correct distance
        print '%s - %s : distance=%d' % (k1, k2, d)
        assert d == distance(n2, n1) # make sure is symmetrical
