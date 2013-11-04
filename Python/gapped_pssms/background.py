#
# Copyright John Reid 2008
#

"""
Code to help use a background model with gapped PSSMs.
"""

import numpy


def _log_c_accumulated(c):
    log_c = numpy.log(c)
    log_c_accumulated = numpy.empty_like(log_c)
    log_c_accumulated[0] = log_c[0]
    log_c_rev_accum = numpy.empty_like(log_c)
    log_c_rev_accum[-1] = log_c[-1]
    for i in xrange(1, len(log_c)):
        log_c_accumulated[i] = log_c_accumulated[i-1] + log_c[i]
        log_c_rev_accum[-i-1] = log_c_rev_accum[-i] + log_c[-i-1]
    return log_c_accumulated, log_c_rev_accum

def forward_backward_log_likelihoods(alpha, beta, c):
    """
    Takes the output of hmm.Model.forward(sequence) (and backward)
    and calculates the log likelihood of all the prefixes
    and suffixes of the sequence.
    """
    log_c_accumulated, log_c_rev_accum = _log_c_accumulated(c)
    log_alpha_sum = numpy.log(alpha.sum(axis=1))
    LL_prefixes = log_alpha_sum - log_c_accumulated
    LL_suffices = numpy.log((beta * alpha).sum(axis=1)) - log_alpha_sum - log_c_rev_accum
    #from IPython.Debugger import Pdb; Pdb().set_trace()
    return LL_prefixes, LL_suffices

def k_mer_log_likelihoods(K, LL, alpha, beta, c):
    """
    Take the output of hmm.Model.forward/backward and use it to calculate the log likelihood of each K-mer.
    """
    LL_prefixes, LL_suffices = forward_backward_log_likelihoods(alpha, beta, c)
    num_k_mers = len(c) - K + 1
    result = numpy.empty((max(0, num_k_mers),))
    if num_k_mers > 0:
        result[0] = LL - LL_suffices[K-1]
        for i in xrange(1, num_k_mers):
            result[i] = LL - LL_prefixes[i-1] - LL_suffices[i+K-1]
    return result
