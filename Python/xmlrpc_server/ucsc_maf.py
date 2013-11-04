#
# Copyright John Reid 2010
#

"""
Interface to UCSC MAF data.
"""

import biopsy.data.ucsc_maf, bx.align.maf, logging

_index = None

from cookbook.cache_decorator import cached_method

@cached_method
def _get_index(genome):
    """
    @return: The index into the MAFs.

    Make sure the index is loaded.

    """
    # maf_files = (biopsy.data.ucsc_maf.maf_file('chr10'),) # replace this with biopsy.data.ucsc_maf.all_maf()
    maf_files = biopsy.data.ucsc_maf.all_maf_files(genome)
    logging.info('Loading index for: %s', ', '.join(maf_files))
    return bx.align.maf.MAFMultiIndexedAccess(
        maf_files,
        keep_open=True,
        parse_e_rows=True,
        use_cache=True
    )


def maf_extract_range_indexed(genome, src, start, end):
    "@return: Get the MAF indexed by the arguments as a sequence of strings."
    blocks = _get_index(genome).get(src, start, end)
    return map(str, blocks)
