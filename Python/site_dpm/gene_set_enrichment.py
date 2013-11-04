#
# Copyright John Reid 2008
#

"""
Code to test for enrichment of gene sets (factors/targets of transcriptional programs)
w.r.t. other gene sets (typically genes in pathways)
"""


import os, biopsy, logging
from itertools import chain, ifilter, imap
from kegg import get_kegg_pathway_info

def read_gene_set_from_file(file):
    "Read a gene set from a file. Each gene is on separate line."
    return set(l.strip() for l in open(file))


def gene_sets_from_dir(dirname):
    "Open each file in the directory and read each line as a gene. Yield (name, gene set) tuples."
    for file in os.listdir(dirname):
        yield file, read_gene_set_from_file(os.path.join(dirname, file))

def kegg_and_inoh_gene_sets():
    "Yield all the (name, gene set) tuples from the KEGG and INOH data."
    return chain(
      gene_sets_from_dir(os.path.join(biopsy.get_data_dir(), 'KEGG', 'pathways')),
      gene_sets_from_dir(os.path.join(biopsy.get_data_dir(), 'INOH', 'pathways')),
    )


class SignificantPValue(object):
    def __init__(self, threshold=.05):
        self.threshold = threshold
        "The threshold at which p-values are called significant"

    def __call__(self, x):
        name, (white_drawn, white, black, draws, p_value)  = x
        return p_value < self.threshold


def test_enrichment(test_set, gene_set, universe):
    """
    Tests the test_set for enrichment relative to the gene_set using the hypergeometric test.
    Returns the probability of seeing at least this many draws in the test set that come from the gene set.

    @arg test_set: The gene set we wish to test for enrichment.
    @arg gene_set: The background set of interest, e.g. genes in a GO category.
    @arg universe: The genes we have selected our gene set from.
    """
    from rpy2.robjects import r
    assert test_set <= universe
    assert gene_set <= universe
    white_drawn = len(test_set.intersection(gene_set))
    white = len(gene_set)
    black = len(universe) - len(gene_set)
    draws = len(test_set)
    #import IPython; IPython.Debugger.Pdb().set_trace()
    kwargs = {'lower.tail' : False}
    logging.debug('%4d black; %4d white; %4d draws; %4d white draws', black, white, draws, white_drawn)
    return white_drawn, white, black, draws, r.phyper(white_drawn-1, white, black, draws, **kwargs)[0]


def test_enrichment_multiple(test_set, gene_sets, universe):
    """
    Expects (name, gene set) tuples from the iterable gene_sets. Calls test_enrichment on each one.
    """
    for name, gene_set in gene_sets:
        yield name, test_enrichment(test_set, universe.intersection(gene_set), universe)



def test_program_target_enrichment(statistics, gene_names, gene_sets, universe):
    """
    Test the targets of the transcriptional programs for enrichment against the gene sets.
    """
    for k in xrange(statistics.num_topics_used):
        test_set = set(gene_names[g] for g in statistics.documents_for_topic(k))
        p_values = [test_enrichment(test_set, gene_set, universe) for gene_set in gene_sets]
        print p_values


class TpEnrichmentTester(object):
    """
    Tests the transcriptional program for enrichment against the factors and the targets, either of which can be None
    """

    def __init__(self, factor_universe, target_universe, p_value_threshold=.05):
        self.factor_universe = set(imap(str, factor_universe))
        self.target_universe = set(imap(str, target_universe))
        self.p_value_threshold = p_value_threshold

    def test_tp_for_factors(self, tp, test_factors):
        "@return: test_drawn, test_size, test_complement_size, draws, p_value"
        return test_enrichment(
            test_factors,
            set(tp.factors),
            self.factor_universe
        )

    def test_tp_for_targets(self, tp, test_targets):
        "@return: test_drawn, test_size, test_complement_size, draws, p_value"
        return test_enrichment(
            test_targets,
            set(tp.targets),
            self.target_universe
        )

    def test_transcriptional_program_factors(self, transcriptional_programs, test_factors):
        test_factors = self.factor_universe.intersection(test_factors)
        for tp in transcriptional_programs:
            test_drawn, test_size, test_complement_size, draws, p_value = self.test_tp_for_factors(tp, test_factors)
            if p_value < self.p_value_threshold:
                yield tp, (test_drawn, test_size, test_complement_size, draws, p_value)

    def test_transcriptional_program_targets(self, transcriptional_programs, test_targets):
        test_targets = self.target_universe.intersection(test_targets)
        for tp in transcriptional_programs:
            test_drawn, test_size, test_complement_size, draws, p_value = self.test_tp_for_targets(tp, test_targets)
            if p_value < self.p_value_threshold:
                yield tp, (test_drawn, test_size, test_complement_size, draws, p_value)


def find_factors_in_enrichment_order(dpm, factors, test_factors, k):
    """
    Looks at the relative enrichment of factors in the posterior of the DPM for each program and yields the indices where the
    test factors come in this ordering.

    This is useful to check if our posterior thresholding is missing some information.
    """
    phi = dpm.exp_phi()
    Phi = dpm.exp_Phi()
    phi_ratio = phi/Phi
    ranking = phi_ratio[k].argsort()[::-1]
    sorted =  [factors[i] for i in ranking]
    for f in test_factors:
        yield sorted.index(f)



if '__main__' == __name__:
    #for name, gene_set in kegg_and_inoh_gene_sets():
    #  print name, len(gene_set)

    universe = set(range(149))
    gene_set = set(range(3))
    test_set = set(range(3))
    print test_enrichment(test_set, gene_set, universe)
