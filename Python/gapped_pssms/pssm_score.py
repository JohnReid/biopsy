#
# Copyright John Reid 2008
#

"""
Code to rank PSSMs by "interesting-ness".

Information content.
Low-order predictability.
Number of sequences with sites.
"""


import numpy

def calculate_first_order_dist(emissions):
    N, M = emissions.shape
    dist = numpy.zeros((M, M))
    for i in xrange(N-1):
        dist += numpy.outer(emissions[i], emissions[i+1])
    dist /= dist.sum()
    return dist

def entropy(dist):
    tmp = -dist * numpy.log(dist)
    tmp[numpy.where(0.0 == dist)] = 0.0
    return tmp.sum()

def calculate_first_order_entropy_score(emissions):
    """
    @return: Score in [0,1] that is ratio of first order entropy to maximum achievable.
    """
    M = emissions.shape[1]
    return entropy(calculate_first_order_dist(emissions))/(2.*numpy.log(M))

def information_contents(emissions):
    M = emissions.shape[1]
    result = emissions * (numpy.log(emissions) + numpy.log(M))
    result[numpy.where(emissions == 0.)] = 0.
    return result.sum(axis=1)

def calculate_information_content_score(emissions):
    """
    @return: Score in [0,1] that is ratio of information content to maximum achievable.
    """
    N, M = emissions.shape
    information_content = information_contents(emissions)
    assert N == len(information_content)
    max_information_content = numpy.log(M)
    return information_content.sum() / max_information_content / N


def geometric_mean(values):
    """
    @return: The geometric mean of the values.
    """
    return numpy.exp(numpy.log(values).sum() / len(values))


def weighted_geometric_mean(values, weights):
    """
    @return: The weighted geometric mean of the values. I.e. the logs of the values are weighted and summed.
    """
    return numpy.exp((numpy.log(values) * weights).sum() / len(values))


def _calculate_emissions(model):
    emissions = numpy.zeros((model.N, model.M))
    for i in xrange(model.N):
        assert model.emissions[i][0] == i
        emissions[i] = model.emissions[i][1]
    return emissions


def show_pssm_score(pssm_filename):
    from parse_gapped_format import parse_models
    models = parse_models(open(pssm_filename))
    for model in models:
        emissions = _calculate_emissions(model)
        first_order_entropy_score = calculate_first_order_entropy_score(emissions)
        first_order_entropy_score **= 2
        information_content_score = calculate_information_content_score(emissions)
        overall_score = geometric_mean((first_order_entropy_score, information_content_score))
        print '%6g %6g %6g' % (first_order_entropy_score, information_content_score, overall_score)

if '__main__' == __name__:
    import sys
    for pssm_filename in sys.argv[1:]:
        print pssm_filename
        show_pssm_score(pssm_filename)
