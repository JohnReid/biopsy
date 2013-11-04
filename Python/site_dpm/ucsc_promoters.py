#
# Copyright John Reid 2009
#

"""
Code to parse the UCSC promoter files Lorenz has provided.
"""


from shared import *
import os, analysis, corebio.seq_io.fasta_io, logging
from corebio.seq import dna_alphabet
from cookbook.dicts import DictOf

ucsc_data_dir = '/home/reid/Data/UCSC/mouse-promoters'
"The directory where the mouse promoters are stored."

fasta_files = [
    os.path.join(ucsc_data_dir, 'promoters_mm9_1000_up_tss.fasta'),
    os.path.join(ucsc_data_dir, 'promoters_mm9_1000_up_tss_2.fasta'),
]
"The FASTA files the unmasked promoters are stored in."


masked_fasta_files = [
    os.path.join(ucsc_data_dir, 'promoters_mm9_1000_up_tss_masked.fasta'),
    #os.path.join(ucsc_data_dir, 'promoters_mm9_1000_up_tss_masked_small.fasta'),
]
"The FASTA files the masked promoters are stored in."


def num_bases_in_promoters(promoters):
    return sum(sum(imap(len, p)) for p in promoters.values())

def get_fasta_files(masked):
    "@return: The FASTA files the sequences are in."
    logging.info('Using %s UCSC promoters' % (masked and 'masked' or 'unmasked'))
    if masked:
        return masked_fasta_files
    else:
        return fasta_files


def get_promoters():
    "Load the promoters into a dict keyed by ensembl gene."
    logging.info('Loading UCSC promoters')
    ucsc_promoters = DictOf(list)
    for fasta_file in get_fasta_files(options.ucsc_use_masked_seqs):
        for seq in corebio.seq_io.fasta_io.iterseq(open(fasta_file), dna_alphabet):
            ucsc, ensembl = seq.description.split()
            ucsc_promoters[ensembl].append(seq)
    logging.info('Loaded %d promoters', len(ucsc_promoters))
    return ucsc_promoters


def do_analysis(ucsc_promoters):
    "Analyse the promoters."
    logging.info('Analysing %d UCSC promoters with %d bases', len(ucsc_promoters), num_bases_in_promoters(ucsc_promoters))
    sequence_analyser = analysis.get_sequence_analyser()
    ucsc_analysis = DictOf(list)
    for ensembl, promoters in ucsc_promoters.iteritems():
        for promoter in promoters:
            ucsc_analysis[ensembl].append(sequence_analyser(promoter))
    return ucsc_analysis


@log_exceptions()
@caching_decorator('ucsc-analysis')
def get_analysis():
    return do_analysis(get_promoters()).copy()


if '__main__' == __name__:
    get_analysis()
