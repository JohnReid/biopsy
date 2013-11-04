#
# Copyright John Reid 2009
#

"""
Code to model transcriptional programs.
"""

from shared import *
import os, csv


def yield_genes_from_csv(f):
    """
    Load the gene names from the file object.
    """
    for l in csv.reader(f,  delimiter=','):
        yield l[0]


class BasicTranscriptionalProgram(object):
    """
    Holds basic information about a transcriptional program.
    """

    def __init__(self, k, factors, targets):
        "Construct a transcriptional program."

        self.k = k
        "The index of this transcriptional program."

        self.factors = factors
        "The factors in this program."

        self.targets = targets
        "The targets of this program."


    def write_files(self, ensembl_names):
        """
        Write the factors and the targets to csv files with their names.
        """
        write_gene_set_with_names(
            open(os.path.join(get_programs_dir(), '%03d-factors.csv' % self.k), 'w'),
            self.factors,
            ensembl_names
        )
        write_gene_set_with_names(
            open(os.path.join(get_programs_dir(), '%03d-targets.csv' % self.k), 'w'),
            self.targets,
            ensembl_names
        )




def tp_from_directory(directory, k):
    """
    Construct a transcriptional program from the factor and target files in the given directory.
    """
    return BasicTranscriptionalProgram(
        k,
        set(yield_genes_from_csv(open(os.path.join(directory, '%03d-factors.csv' % k)))),
        set(yield_genes_from_csv(open(os.path.join(directory, '%03d-targets.csv' % k))))
    )


def tp_from_dpm_summary(summariser, factor_universe, target_universe, k):
    """
    Construct a transcriptional program from the summary of a DPM.
    """
    return BasicTranscriptionalProgram(
        k,
        set(str(factor_universe[w]) for w in summariser.statistics.words_for_topic(k)),
        set(str(target_universe[d]) for d in summariser.statistics.documents_for_topic(k))
    )
