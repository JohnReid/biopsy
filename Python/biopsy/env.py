#
# Copyright John Reid 2006
#

import os
import biopsy

def get_data_dir():
    return biopsy.Environment.data_dir

def get_biopsy_data_dir():
    return os.path.join(get_data_dir(), 'biopsy')

def get_aliases_dir():
    return os.path.join(get_data_dir(), 'Aliases')
