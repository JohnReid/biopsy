#
# Copyright John Reid 2010
#

"""
Package to handle various data sources stored in the biopsy data directory.
"""

import os

_potential_data_dirs = (
    '/home/john/Data',
    '/home/reid/Data',
)


def data_dir():
    "@return: The root data directory."
    for dir in _potential_data_dirs:
        if os.path.exists(dir):
            return dir
    raise RuntimeError('Do not know where data directory is on this box.')
