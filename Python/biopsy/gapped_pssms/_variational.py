#
# Copyright John Reid 2006
#

import numpy, numpy.random, scipy.special, math
from _maths import *

class VariationalModel( object ):
    """A variational model over hidden variables in gapped pssm problem
    """

    only_sites = True

    def __init__(
            self,
            seqs,
            K,
            alpha,
            phi,
            varphi,
    ):
        """Initialise a variational model

        seqs: the sequences
        K: the size of the PSSM
        """
        self.seqs = seqs
        self.N = len( self.seqs )
        self.alpha = alpha
        self.phi = phi
        self.K = K
        self.varphi = varphi
        self.init_var_params()

    def init_var_params( self ):
        """Initialise variational parameters"""
        # parameter for theta's
        self.omega = numpy.empty( (self.K+1, 4), numpy.float64 )
        for j in xrange( self.K + 1 ):
            if 0 == j: self.omega[j,:] = self.varphi
            else: self.omega[j,:] = sample_from_dirichlet( self.phi )
        # parameter for gamma
        self._lambda = numpy.ones( 2, dtype = numpy.float64 )
        # parameter for h
        self.eta = sample_from_dirichlet( [ 1.0 ] * (self.K-1) )
        # parameter for g's
        self.mu = sample_from_dirichlet( [ 1.0 ] * self.N )
        # parameter for s's
        self.nu = [
                sample_from_dirichlet( [ 1.0 ] * ( len(s) - self.K ) )
                for s in self.seqs
        ]
        self._check()

        # update cached values
        self._recalc_after_omega_changes()

    def expected_pssm( self ):
        """The expected pssm distribution - including base distribution as index 0"""
        e_pssm = numpy.array( self.omega, copy = True )
        for row in e_pssm:
            row /= numpy.sum( row )
        return e_pssm

    def pssm_entropy_per_base( self ):
        """The entropy of the theta distribution over the pssm bases"""
        return sum(
                [
                        expected_entropy_under_dirichlet( self.omega[j + 1,:] )
                        for j in xrange( self.K )
                ]
        ) / self.K

    def bg_entropy( self ):
        """The entropy of the background distribution over the background bases"""
        return expected_entropy_under_dirichlet( self.omega[0,:] )

    def start_entropy_per_sequence( self ):
        """The entropy of the start distribution"""
        return sum( [ discrete_entropy( nu ) for nu in self.nu ] ) / self.N

    def has_gap_entropy_per_sequence( self ):
        """The entropy of the has gap distribution"""
        return numpy.sum( [ binary_entropy( mu ) for mu in self.mu ] ) / self.N

    def position_entropy( self ):
        """The entropy of the position distribution"""
        return discrete_entropy( self.eta )

    def _check( self ):
        """Checks that everything is ship shape"""
        for n in xrange( self.N ):
            assert finite_and_a_number( self.mu[n] )
            assert (self.mu[n] > 0).all()
            assert finite_and_a_number( self.nu[n] )
            assert (self.nu[n] > 0).all()
        assert finite_and_a_number( self.eta )
        assert (self.eta > 0).all()
        assert finite_and_a_number( self.omega )
        assert (self.omega > 0).all()
        assert self.K > 0
        assert self.eta.shape == (self.K-1,)

    def __str__( self ):
        return 'N: %d\nK: %d\nalpha: %s\nphi: %s\nvarphi: %s\nomega:\n%s\n' \
                'lambda: %s\neta: %s\nmu: %s' % (
                self.N,
                self.K,
                str( self.alpha ),
                str( self.phi ),
                str( self.varphi ),
                str( self.omega ),
                str( self._lambda ),
                str( self.eta ),
                str( self.mu ),
        )

    def _expectation_log_gamma( self ):
        return (
                scipy.special.digamma( self._lambda[0] )
                - scipy.special.digamma( numpy.sum( self._lambda ) )
        )

    def _expectation_log_1_minus_gamma( self ):
        return (
                scipy.special.digamma( self._lambda[1] )
                - scipy.special.digamma( numpy.sum( self._lambda ) )
        )

    def _recalc_after_omega_changes( self ):
        # print 'omega\n%s' % str(self.omega)
        self.log_p_x_given_r = numpy.empty_like( self.omega )
        self.p_x_given_r = numpy.empty_like( self.omega )
        for r in xrange( self.K + 1 ):
            for x in xrange( 4 ):
                t = numpy.sum(self.omega[r,:])
                self.log_p_x_given_r[r,x] = (
                        scipy.special.digamma( self.omega[r,x] )
                        - scipy.special.digamma( t )
                )
                self.p_x_given_r[r,x] = self.omega[r,x] / t

    class Combinations( object ):
        """Generator of all combinations of i, s, g, h with their likelihoods
        under the variational distribution

        Attributes:
                s: start position of binding site
                q_s: likelihood of this start position
                h: gap position
                q_h: likelihood of this gap position
                g: binding site has gap?
                q_g: likelihood of g
                i: index in sequence
                x: base
                r: index into pssm
        """
        def __init__( self, model ):
            self.model = model

        def __call__( self, n, only_sites = False ):
            """A generator for all combinations

            n: the sequence for the combinations to be generated over
            only_sites: if true then positions where r = 0 are not returned
            """
            model = self.model
            g = [ 1.0 - model.mu[n], model.mu[n] ]
            for self.s, self.q_s in enumerate( model.nu[n] ):
                if only_sites: i_generator = xrange( self.s, self.s + model.K + 1 )
                else: i_generator = xrange( len( model.seqs[n] ) )
                for self.i in i_generator:
                    self.x = model.seqs[n][self.i]

                    j = self.i - self.s
                    outside_pssm = j < 0 or j > model.K

                    for self.h, self.q_h in enumerate( model.eta ):
                        for self.g, self.q_g in enumerate( g ):
                            # gap position
                            if self.g: gap_position = self.h
                            else: gap_position = model.K - 1

                            # calculate r - index into pssm - 0 is background
                            if outside_pssm: self.r = 0
                            elif j == gap_position + 1: self.r = 0
                            elif j > gap_position: self.r = j
                            else: self.r = j + 1
                            assert 0 <= self.r < model.K + 1

                            if only_sites and self.r == 0: continue

                            yield self

    def log_likelihood( self ):
        LL = 0.0
        for n in xrange( self.N ): # for each sequence
            likelihoods = numpy.zeros( len( self.seqs[n] ), dtype = numpy.float64 )
            for c in self.Combinations( self )( n, only_sites = False ):
                likelihoods[c.i] += (
                        c.q_s * c.q_h * c.q_g
                        * self.p_x_given_r[c.r,c.x]
                )
            LL += numpy.sum( numpy.log( likelihoods ) )
        return LL

    def update( self ):
        """One round of variational updates"""
        self._update_omega()
        self._update_mu_nu_eta()
        self._update_lambda()

    def _update_mu_nu_eta( self ):
        """Update mu, nu and eta, the variational parameters for g, s and h"""
        log_h = numpy.zeros_like( self.eta )
        for n in xrange( self.N ):
            log_g = numpy.array(
                    [
                            self._expectation_log_gamma(),
                            self._expectation_log_1_minus_gamma(),
                    ],
                    dtype = numpy.float64
            )
            log_s = numpy.zeros_like( self.nu[n] )
            for c in self.Combinations( self )( n, only_sites = self.only_sites ):
                g_weight = c.q_s * c.q_h # weight in expectation for g
                s_weight = c.q_g * c.q_h # weight in expectation for s
                h_weight = c.q_g * c.q_s # weight in expectation for h
                expectation = self.log_p_x_given_r[c.r,c.x]
                log_g[c.g] += g_weight * expectation
                log_s[c.s] += s_weight * expectation
                log_h[c.h] += h_weight * expectation

            # update mu
            self.mu[n] = probabilities_from_logs( log_g )[1]

            # update nu
            self.nu[n] = probabilities_from_logs( log_s )

        # update eta
        self.eta = probabilities_from_logs( log_h )
        self._check()

    def _update_lambda( self ):
        """Update lambda, the variational parameters for g"""
        t = numpy.sum( self.mu )
        self._lambda[0] = self.alpha[0] + t
        self._lambda[1] = self.alpha[1] + self.N - t

    def _update_omega( self ):
        """Update omega, the variational parameters for pi and theta"""
        log_theta = numpy.empty_like( self.omega )
        for j in xrange( self.K + 1 ):
            if 0 == j: prior = self.varphi
            else: prior = self.phi
            for x in xrange( 4 ):
                log_theta[j,x] = prior[x]
        for n in xrange( self.N ):
            for c in self.Combinations( self )( n, only_sites = self.only_sites ):
                log_theta[c.r,c.x] += c.q_s * c.q_h * c.q_g

        # update omega
        #for j in xrange( self.K + 1 ):
        #       for x in xrange( 4 ):
        #               log_theta[r,x] *= self.log_p_x_given_r[r,x]
        self.omega = log_theta

        # update cached values
        self._recalc_after_omega_changes()
        self._check()

    def most_likely( self ):
        return (
                [
                        numpy.argmax( nu )
                        for nu in self.nu
                ],
                [
                        mu > 0.5
                        for mu in self.mu
                ],
                numpy.argmax( self.eta )
        )

    class Learner( object ):
        def __init__(
                self,
                model,
                max_updates = None,
                convergence_epsilon = 1e-6
        ):
            self.model = model
            self.LL = [ self.model_LL() ]
            self.max_updates = max_updates
            self.convergence_epsilon = convergence_epsilon
            self.gain = None

        def model_LL( self ): return self.model.log_likelihood()
        def update_model( self ): self.model.update()

        def update( self ):
            """Do one update"""
            self.update_model()
            LL = self.model_LL()
            self.gain = (LL - self.LL[-1]) / numpy.fabs(LL)
            self.LL.append( LL )
            return LL

        def _terminate( self ):
            # have we performed the maximum number of updates?
            if self.max_updates and self.max_updates < len( self.LL ):
                return True

            # did the last 2 updates not change the LL significantly
            if (
                    len( self.LL ) > 2
                    and
                    math.fabs( (self.LL[-1] - self.LL[-2]) / self.LL[-2] ) < self.convergence_epsilon
                    and
                    math.fabs( (self.LL[-2] - self.LL[-3]) / self.LL[-3]  ) < self.convergence_epsilon
            ):
                return True

            return False

        def __call__( self ):
            while not self._terminate():
                yield self.update()

        def learn(
                self,
                callback = None
        ):
            for LL in self():
                if callback: callback( self )
