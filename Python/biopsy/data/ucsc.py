#
# Copyright John Reid 2010, 2011, 2012
#

"""
Access UCSC data stored locally.
"""

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import gzip, os, logging

logger = logging.getLogger(__name__)

ucsc_data_dir = "/home/john/Data/UCSC"

def full_path(filename):
    "@return: The full path for the filename."
    return os.path.join(ucsc_data_dir, filename)

def upstream_5000(genome='mm9'):
    return upstream(genome, length=5000)

def upstream(genome='mm9', length=5000):
    "@return: Yield 5000bp sequences upstream of the TSS."
    filename = full_path(os.path.join('goldenPath', genome, 'bigZips', 'upstream%d.fa.gz' % length))
    logger.info('Loading upstream sequences from %s', filename)
    for record in SeqIO.parse(
        gzip.open(filename, 'rb'),
        "fasta",
        generic_dna
    ):
        yield record


def chromosome_filename(genome, chromsome):
    "@return: The filename of the chromosome."
    return full_path(os.path.join('goldenPath', genome, 'chromosomes', '%s.fa.gz' % chromsome))


class Genome(dict):
    """
    @return: A dictionary mapping chromosome names to sequences.
    """
    def __init__(self, genome_name):
        self.genome_name = genome_name
        "The name of the genome, e.g. mm8."
        
    def __missing__(self, chromosome):
        filename = chromosome_filename(self.genome_name, chromosome)
        logger.info('Loading %s %s from %s', self.genome_name, chromosome, filename)
        seqs = list(SeqIO.parse(
            gzip.open(filename),
            "fasta",
            generic_dna)
        )
        if 1 != len(seqs):
            raise RuntimeError('Expecting exactly one sequence in %s' % filename)
        self[chromosome] = seqs[0]
        return seqs[0]

if '__main__' == __name__:
    for record in upstream_5000():
        print record
