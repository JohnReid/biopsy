/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/contingency_homogeneity.h"
#include "bio/raii.h"
#include "bio/gsl.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_pow_int.h>

BIO_NS_START

double
gamma_exponential_prior_integrand(double c, void * void_params)
{
	const GammaParams & params = *static_cast<GammaParams *>(void_params);
	double multiplicative_term = 1;
	for (unsigned j = 1; j <= params.k; ++j)
	{
		multiplicative_term *= (c + (j - 1));
	}
	//work in logarithms
	double ln_result = params.m * gsl_sf_log(c);
	for (unsigned j = 1; j <= params.k; ++j)
	{
		ln_result -= gsl_sf_log(c + (j - 1));
	}
	ln_result += (*params.prior).get_ln_prior(c); //multiply by prior

	//this may underflow so check...
	return ln_result < GSL_LOG_DBL_MIN ? 0 : gsl_sf_exp(ln_result);
}

double
ContingencyClusterPrior::get_ln_prior(double c)
{
	return gsl_sf_log((*this)(c));
}

ExponentialPrior::ExponentialPrior(double mean)
: mean(mean)
{
}

double
ExponentialPrior::operator()(double c)
{
	const double exp_param = -(mean * c);
	return exp_param < GSL_LOG_DBL_MIN ? 0 : (mean * gsl_sf_exp(exp_param));
}

double
ExponentialPrior::get_ln_prior(double c)
{
	return gsl_sf_log(mean) - (mean * c);
}

double
contingency_calculate_gamma(
	unsigned m,
	unsigned k,
	ContingencyClusterPrior * prior)
{
	static const unsigned workspace_size = 1000;
	static RAII<gsl_integration_workspace *> w(gsl_integration_workspace_alloc(workspace_size));

	double result, error;

	GammaParams params;
	params.m = m;
	params.k = k;
	params.prior = prior;

	gsl_function F;
	F.function = &gamma_exponential_prior_integrand;
	F.params = &params;

	gsl_integration_qagiu(
		&F, //function
		0, //a
		0, //epsabs
		1e-7, //epsrel
		workspace_size, //limit
		w, //workspace
		&result, //result
		&error);  //error

	//std::cout << "result          = " << result << std::endl;
	//std::cout << "estimated error = " << error << std::endl;
	//std::cout << "intervals       = " << w->size << std::endl;

	return result;
}

BIO_NS_END
