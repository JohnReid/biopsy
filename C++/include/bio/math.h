

#ifndef BIO_MATH_H_
#define BIO_MATH_H_

#include "bio/defs.h"
#include "bio/gsl.h"

#include <boost/io/ios_state.hpp>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include <numeric>


BIO_NS_START

/** Calculates the log of n choose the observed counts.  */
template < typename ObsIt >
double
calc_ln_n_choose(
	ObsIt obs_begin,
	ObsIt obs_end)
{
	double result = 0.0;

	unsigned n = 0;
	while (obs_end != obs_begin)
	{
		result -= gsl_sf_lnfact(*obs_begin);
		n += *obs_begin;

		++obs_begin;
	}

	result += gsl_sf_lnfact(n);

	return result;
}


/** Calculates a factor involved in the likelihood of a multinomial distribution with a dirichlet prior. */
template < typename ObsIt, typename AlphaIt >
double
calc_ln_gamma_factor(
	ObsIt obs_begin,
	ObsIt obs_end,
	AlphaIt alpha_begin)
{
	double ln_pi_gamma_alpha_obs = 0.0; //work in logarithms
	double ln_pi_gamma_alpha = 0.0; //work in logarithms
	double sum_alpha = 0.0;
	double sum_obs = 0.0;
	while (obs_end != obs_begin)
	{
		sum_alpha += *alpha_begin;
		sum_obs += *obs_begin;
		ln_pi_gamma_alpha += gsl_sf_lngamma(*alpha_begin);
		ln_pi_gamma_alpha_obs += gsl_sf_lngamma(*alpha_begin + *obs_begin);

		++obs_begin;
		++alpha_begin;
	}

	double ln_result =
		gsl_sf_lngamma(sum_alpha)
		- ln_pi_gamma_alpha
		+ ln_pi_gamma_alpha_obs
		- gsl_sf_lngamma(sum_alpha + sum_obs);

	return ln_result;
}

/** Calculates the log of the likelihood of seeing the observed counts under a dirichlet prior. */
template < typename ObsIt, typename AlphaIt >
double
calc_multinomial_ln_likelihood_dirichlet_prior(
	ObsIt obs_begin,
	ObsIt obs_end,
	AlphaIt alpha_begin)
{
	return calc_ln_gamma_factor(obs_begin, obs_end, alpha_begin) + calc_ln_n_choose(obs_begin, obs_end);
}


/** Calculates the log of the likelihood of seeing the observed counts under a uniform multinomial distribution. */
template < typename ObsIt >
double
calc_multinomial_ln_likelihood_uniform_dist(
	ObsIt obs_begin,
	ObsIt obs_end)
{
	const unsigned n = std::accumulate(obs_begin, obs_end, unsigned(0));
	if (0 == n)
	{
		return 0.0;
	}

	const double ln_n_choose = calc_ln_n_choose(obs_begin, obs_end);
	return ln_n_choose - n * gsl_sf_log(double(obs_end - obs_begin));
}


/** Calculates the log of the likelihood of seeing the observed counts under a multinomial distribution with given p's. */
template < typename ObsIt, typename PIt >
double
calc_multinomial_ln_likelihood(
	ObsIt obs_begin,
	ObsIt obs_end,
	PIt p_begin)
{
	const unsigned n = std::accumulate(obs_begin, obs_end, unsigned(0));
	if (0 == n)
	{
		return 0.0;
	}

	double result = calc_ln_n_choose(obs_begin, obs_end);
	while (obs_end != obs_begin)
	{
		result += *obs_begin * gsl_sf_log(*p_begin);

		++obs_begin;
		++p_begin;
	}

	return result;
}


/** Calculates the ratio of evidence in favour of a uniform distribution over one with a dirichlet prior
defined by the alphas. */
template < typename ObsIt, typename AlphaIt >
double
calc_ln_multinomial_uniform_evidence(
	ObsIt obs_begin,
	ObsIt obs_end,
	AlphaIt alpha_begin)
{
	const double ln_likelihood_uniform_dist = calc_multinomial_ln_likelihood_uniform_dist(obs_begin, obs_end);
	const double ln_likelihood_dirichlet_prior = calc_multinomial_ln_likelihood_dirichlet_prior(obs_begin, obs_end, alpha_begin);
	return ln_likelihood_uniform_dist - ln_likelihood_dirichlet_prior;
}

/** Calculates the ratio of evidence in favour of a particular multinomial distribution (defined by the p's)
over one with a dirichlet prior defined by the alpha's. */
template < typename ObsIt, typename AlphaIt, typename PIt >
double
calc_ln_multinomial_evidence(
	ObsIt obs_begin,
	ObsIt obs_end,
	AlphaIt alpha_begin,
	PIt p_begin)
{
	return
		calc_multinomial_ln_likelihood(obs_begin, obs_end, p_begin)
		- calc_multinomial_ln_likelihood_dirichlet_prior(obs_begin, obs_end, alpha_begin);
}


BIO_NS_END


#endif //BIO_MATH_H_

