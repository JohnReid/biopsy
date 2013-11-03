/** @file
 * Based on the paper:
 * Quintana, F.A. (1998) ``Nonparametric Bayesian analysis for assessing homogeneity in $ k\times l$
 * contingency tables with fixed right margin totals''.
 * Also as Technical Report PUC/FM-96/7. Journal of the American Statistical Association, 93(443), 1140-1149.
 * Available as compressed postscript.
 * http://www.mat.puc.cl/~quintana/hcttr.ps.gz
 */

#ifndef BIO_CONTINGENCY_HOMOGENEITY_H_
#define BIO_CONTINGENCY_HOMOGENEITY_H_

#include "bio/defs.h"
#include "bio/partition.h"

#include <boost/assign/list_of.hpp>

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>

#include <numeric>
#include <stdexcept>

BIO_NS_START

typedef std::vector<double> gamma_vec_t;
extern gamma_vec_t contingency_gamma_exponential_prior;

template <class LambdaIt>
double
calculate_ln_D(
	LambdaIt lambda_begin,
	LambdaIt lambda_end)
{
	double result = gsl_sf_lngamma(std::accumulate(lambda_begin, lambda_end, double(0)));
	for (LambdaIt i = lambda_begin; lambda_end != i; ++i)
	{
		result -= gsl_sf_lngamma(*i);
	}

	return result;
}

/** Calculate natural logarithm of a term for one partition in Lm(n). */
template <class PartitionIt, class CountIt, class LambdaIt>
double
contingency_calculate_ln_L_term(
	PartitionIt part_it,
	unsigned m,
	CountIt n_begin,
	CountIt n_end,
	LambdaIt lambda_begin,
	LambdaIt lambda_end)
{
	double ln_result = 1.0;

#ifdef BIO_CONTINGENCY_PRINT_PARTITIONS
	//can only do this for smallish numbers
	if (26 >= part_it.M.size())
	{
		std::vector<char> v;
		for (unsigned i = 0; i != part_it.M.size(); ++i)
		{
			v.push_back('a' + i);
		}
		std::cout << *part_it[v] << '\n';
	}
#endif //BIO_CONTINGENCY_PRINT_PARTITIONS

	//for each set in the partition
	for (unsigned i = 0; i < m; ++i)
	{
		//calculate lambda + n
		typedef std::vector<double> lambda_vec_t;
		lambda_vec_t lambda_plus_n;
		{
			std::copy(lambda_begin, lambda_end, std::inserter(lambda_plus_n, lambda_plus_n.begin()));
			unsigned e_i = 0;
			unsigned j = 0;
			//for each index
			for (CountIt n = n_begin; n_end != n; ++j, ++n)
			{
				//is this j in this partition?
				if (i == part_it.M[j])
				{
					for (unsigned k = 0; k != (*n).size(); ++k)
					{
						lambda_plus_n[k] += (*n)[k];
					}
					if (1 < e_i)
					{
						ln_result += gsl_sf_log(double(e_i)); //the factorial term
					}
					++e_i;
				}
			}
		}

		ln_result += calculate_ln_D(lambda_begin, lambda_end);
		ln_result -= calculate_ln_D(lambda_plus_n.begin(), lambda_plus_n.end());
	}

	return ln_result;
}

/** Calculate Lm(n). */
template <class CountIt, class LambdaIt>
double
contingency_calculate_L(
	unsigned m,
	CountIt n_begin,
	CountIt n_end,
	LambdaIt lambda_begin,
	LambdaIt lambda_end)
{
	unsigned set_size = unsigned(n_end - n_begin);
	BOOST_ASSERT(unsigned(lambda_end - lambda_begin) == set_size);

	double result = 0;

	//iterate over all the partitions
	partition::iterator_k part_it(set_size, m);
	try
	{
		while (true)
		{
			const double ln_L_term =
				contingency_calculate_ln_L_term(
					part_it,
					m,
					n_begin,
					n_end,
					lambda_begin,
					lambda_end);

			//check for underflow - in which case it will be 0.
			if (ln_L_term > GSL_LOG_DBL_MIN)
			{
				result += gsl_sf_exp(ln_L_term);
			}

			++part_it;
		}
	}
	catch (std::overflow_error&) //how the partition code determines there are no more partitions!
	{
	}

	return result;
}

template <class Value>
Value
factorial(Value value)
{
	Value result = Value(1);
	while (value != 0)
	{
		result *= value;
		--value;
	}
	return result;
}

typedef double(* contingency_cluster_prior)(double c, void * void_params);

struct ContingencyClusterPrior
{
	virtual double operator()(double c) = 0;
	virtual double get_ln_prior(double c);
};

struct ExponentialPrior : ContingencyClusterPrior
{
	double mean;
	ExponentialPrior(double mean = 1.0);
	double operator()(double c);
	double get_ln_prior(double c);
};

struct GammaParams
{
	unsigned m;
	unsigned k;
	ContingencyClusterPrior * prior;
};

double
gamma_exponential_prior_integrand(double c, void * void_params);

double
contingency_calculate_gamma(
	unsigned m,
	unsigned k,
	ContingencyClusterPrior * prior);


template <class GammaIt, class CountIt, class LambdaIt>
double
contingency_calculate_bayes_factor(
	GammaIt gamma_begin,
	GammaIt gamma_end,
	CountIt n_begin,
	CountIt n_end,
	LambdaIt lambda_begin,
	LambdaIt lambda_end)
{
	const unsigned k = unsigned(n_end == n_begin ? 0 : n_begin->end() - n_begin->begin());

	double result = contingency_calculate_L(1, n_begin, n_end, lambda_begin, lambda_end);

	result *= (1 - factorial(k-1) * (*gamma_begin));

	result /= factorial(k-1);

	double sum = 0;
	for (unsigned i = 2; k + 1 != i; ++i)
	{
		sum += *(gamma_begin + i - 1) * contingency_calculate_L(i, n_begin, n_end, lambda_begin, lambda_end);
	}
	result /= sum;

	return result;
}

/** Contains the data for a contingency table and defines types. */
template <unsigned K, unsigned L>
struct ContingencyTable
{
	typedef boost::array<unsigned, K> row_t;
	typedef boost::array<row_t, L> table_t;
	typedef boost::array<double, L> lambda_vec_t;

	table_t data;

	static unsigned get_k() { return K; }
	static unsigned get_l() { return L; }
};

template <unsigned K, unsigned L>
std::ostream &
operator<<(std::ostream & os, const ContingencyTable<K, L> & table)
{
	for (unsigned l = 0; L != l; ++l)
	{
		for (unsigned k = 0; K != k; ++k)
		{
			os << table.data[l][k] << '\t';
		}
		os << '\n';
	}
	return os;
}



BIO_NS_END


#endif //BIO_CONTINGENCY_HOMOGENEITY_H_
