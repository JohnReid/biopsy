#
# Copyright John Reid 2009
#

"""
Code to generate statistics for transcriptional programs.
"""

from itertools import imap, cycle
from infpy.bootstrap import *



def create_chi_squared_statistic(gene_sets, gene_universe):
    """
    @return: A function that takes a set of genes and sums chi squared statistics against the gene sets
    """
    # make everything into sets
    gene_universe = set(gene_universe)
    gene_sets = [set(gs).intersection(gene_universe) for gs in gene_sets]

    # calculate proportions one time
    proportions = [float(len(gs)) / float(len(gene_universe)) for gs in gene_sets]

    # define statistic function
    def chi_squared_stat(test_gene_set):
        """
        Sums chi squared tests for test_gene_set against each gene set in gene_sets.
        """
        if not test_gene_set:
            return 0.

        def chi_squared(gs, p):
            "@return: Chi squared stat for test set against one gene set."
            expected = p * float(len(test_gene_set))
            observed = len(gs.intersection(test_gene_set))
            return ((observed-expected)**2)/expected

        return sum(map(chi_squared, gene_sets, proportions))

    # return it
    return chi_squared_stat





if '__main__' == __name__:
    test_universe = tester.target_universe
    test_sets = [tp.tp_targets for tp in transcriptional_programs]
    chi_squared_stat = create_chi_squared_statistic(highly_expressed_genes, test_universe)
    bootstrap_stats = calculate_bootstrap_statistics(
        generate_bootstrap_samples(1000, test_universe, map(len, test_sets)),
        chi_squared_stat
    )
    for k, ts in enumerate(test_sets):
        stat = chi_squared_stat(ts)
        p_value = bootstrap_p_value(bootstrap_stats, stat)
        if p_value < .05:
            print k, p_value, stat
