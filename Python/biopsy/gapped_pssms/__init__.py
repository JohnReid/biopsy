#
# Copyright John Reid 2006
#

import numpy, numpy.random, scipy.special, math

from _maths import *
from _generate import *
from _generate_2 import *
from _variational import *
from gapped_pssms import *
from weblogo import *


#
# Try to import C++ part of module if installed.
#
try:
    from _gapped_pssms import *
    #
    # The c implementation does not hold the data as numpy arrays
    # so provide some functions to create numpy arrays from the data
    #
    def _gapped_pssm_alpha_array( model ):
        """Beta prior parameters for gamma: the likelihood of a gap"""
        return numpy.array(
                [
                        model.alpha( i )
                        for i in xrange( 2 )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.alpha_array = _gapped_pssm_alpha_array

    def _gapped_pssm_varphi_array( model ):
        "Dirichlet prior parameters for background distribution"
        return numpy.array(
                [
                        model.varphi( i )
                        for i in xrange( 4 )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.varphi_array = _gapped_pssm_varphi_array

    def _gapped_pssm_phi_array( model ):
        "Dirichlet prior parameters for pssm distribution"
        return numpy.array(
                [
                        model.phi( i )
                        for i in xrange( 4 )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.phi_array = _gapped_pssm_phi_array

    def _gapped_pssm_lambda_array( model ):
        "Variational parameter for gamma"
        return numpy.array(
                [
                        model.lambda_( i )
                        for i in xrange( 2 )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.lambda_array = _gapped_pssm_lambda_array

    def _gapped_pssm_eta_array( model ):
        "Variational parameter for location of the gap"
        return numpy.array(
                [
                        model.eta( i )
                        for i in xrange( model.K - 1 )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.eta_array = _gapped_pssm_eta_array

    def _gapped_pssm_mu_array( model ):
        "Variational parameter for g: has_gap variable"
        return numpy.array(
                [
                        model.mu( i )
                        for i in xrange( model.N )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.mu_array = _gapped_pssm_mu_array

    def _gapped_pssm_omega_array( model ):
        "Variational parameters for background and pss distributions"
        return numpy.array(
                [
                        [
                                model.omega( r, x )
                                for x in xrange( 4 )
                        ]
                        for r in xrange( model.K+1 )
                ],
                dtype = numpy.float64
        )
    VariationalModel_C.omega_array = _gapped_pssm_omega_array

    def _gapped_pssm_nu_sequence( model ):
        "Variational parameters for start positions of sites"
        return [
                numpy.array(
                        [
                                model.nu( n, i )
                                for i in xrange( 2 * (model.sequence_length( n ) - model.K) )
                        ],
                        dtype = numpy.float64
                )
                for n in xrange( model.N )
        ]
    VariationalModel_C.nu_sequence = _gapped_pssm_nu_sequence

except ImportError:
    import warnings
    warnings.warn('Could not import C++ gapped PSSM module')
