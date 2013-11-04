#
# Copyright John Reid 2010
#

"""
Use UCSC MAF files.
"""


import ucsc, bx.align.maf, glob, os

golden_path_dir = '/opt/genomics/UCSC/goldenPath'


def full_path(filename):
    "@return: The full path for the filename."
    return os.path.join(golden_path_dir, filename)


_dirs_for_genomes = {
    'mm9'    : full_path('mm9/multiz30way/maf'),
    'hg19'   : full_path('hg19/multiz46way/maf'),
}
"Maps genomes to their alignment directories."


def alignment_dir(genome):
    "@return: The directory that the alignment files for the genome are stored in."
    if genome in _dirs_for_genomes:
        return _dirs_for_genomes[genome]
    else:
        raise RuntimeError('Do not have alignment directory for genome %s' % genome)


def alignment_file(genome, file):
    "@return: The full pathname of the alignment file for the genome."
    return os.path.join(alignment_dir(genome), file)


def all_maf_files(genome):
    "@return: All MAF files."
    return glob.glob(alignment_file(genome, '*.maf.bz2'))


def maf_file(genome, chr):
    "@return: The MAF file for the chromosome."
    return alignment_file(genome, '%s.maf.bz2' % chr)


if '__main__' == __name__:
    genome = 'mm9'
    maf_files = all_maf_files(genome)

    index = bx.align.maf.MultiIndexed(
        (maf_file(genome, 'chr10'),),
        keep_open=True,
        parse_e_rows=True,
        use_cache=True
    )

    src = '%s.chr10' % genome
    start = 10000000
    end = 10000300
    blocks = index.get(src, start, end)

    for block in blocks:
        #print type(block)
        print block
