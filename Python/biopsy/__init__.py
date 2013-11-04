#
# Copyright John Reid 2006, 2013
#

from _biopsy import *
init()

import os

import transfac
from transfac import db, DbRef

from env import *
from binding_hit import *
from chromosome import *
from families import *
from gene_locations import *
from graph import *
from lcs import *
from pair_analysis import *
from paralogs import *
from pssm import *
from remo import *
from remome import *
from sequence import *
from test_case import *
from web import *

import pickle, cPickle

def serialisation_dir():
    import os
    _dir = os.path.join(biopsy.get_data_dir(), 'biopsy/')
    if not os.access(_dir, os.X_OK): os.makedirs(_dir)
    return _dir

def load(
        f,
        display = None
):
    """Loads an object from the file (specified by filename or file object)"""
    if None != display: print 'LOAD; %s; %s' % ( f, display )
    if isinstance( f, str ): f = open( f, 'rb' )
    return cPickle.load( f )

def save(
        obj,
        f,
        display = None,
        protocol = cPickle.HIGHEST_PROTOCOL
):
    """Loads an object from the file (specified by filename or file object)"""
    if None != display: print 'SAVE; %s; %s' % ( f, display )
    if isinstance( f, str ): f = open( f, 'wb' )
    pickle.dump( obj, f, protocol = protocol )
    # try:
    #       cPickle.dump( obj, f, protocol = protocol )
    # except:
    #       pickle.dump( obj, f, protocol = protocol )


_base_url = 'http://www.biobase-international.com/'
_version = '10.1'

#def get_pssm_url( pssm_acc ):
#       return \
#               _base_url \
#               + 'cgi-bin/biobase/transfac/' \
#               + _version \
#               + '/bin/getTFProf.cgi?' \
#               + str( pssm_acc )

def find_transfac_pssms_by_name( name, accessions = None ):
    """Returns a list of those pssms whose name contains the given name"""
    if None == accessions:
        accessions = get_transfac_pssm_accessions(
                get_default_transfac_pssm_filter()
        )
    return [
            p
            for p
            in accessions
            if -1 != get_transfac_pssm_name( p ).find( name )
    ]


def _pssm_numpy_repr(pssm):
    """Return a numpy representation of the frequencies of a PSSM.
    """
    import numpy
    return numpy.array([[col.get_freq(b) for b in xrange(4)] for col in pssm])


def _information_content_matrix(matrix):
    """The information content of the matrix of a PSSM in bits."""
    from numpy import log2
    B = matrix.shape[1]
    return _safe_x_log2_x(matrix).sum() + matrix.shape[0] * log2(B)


def _information_content(pssm):
    """The information content of a PSSM in bits."""
    return _information_content_matrix(pssm.as_numpy())


def _safe_x_log2_x(x):
    """x * log2(x), with zeros where x is zero
    """
    from numpy import log2
    result = log2(x) * x
    result[0==x] = 0.
    return result


# Add to the Pssm class
Pssm.as_numpy = _pssm_numpy_repr
Pssm.information_content = _information_content



