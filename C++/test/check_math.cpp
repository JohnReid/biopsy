

#include "bio/math.h"
#include "bio/iterator.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;

#include <gsl/gsl_sf_exp.h>


//#define VERBOSE_CHECKING




typedef vector< unsigned > obs_vec_t;

const std::vector< double > alpha(4, 1.0);
const vector< obs_vec_t > obs_vec = list_of< vector< unsigned > >
	(list_of< unsigned > (0) (0) (0) (0))
	(list_of< unsigned > (0) (1) (1) (0))
	(list_of< unsigned > (0) (2) (0) (0))
	(list_of< unsigned > (1) (1) (1) (1))
	(list_of< unsigned > (4) (4) (4) (4))
	(list_of< unsigned > (0) (0) (0) (4))
	(list_of< unsigned > (0) (0) (0) (8))
	(list_of< unsigned > (1) (1) (1) (13))
	(list_of< unsigned > (100) (97) (103) (101))
	;

const vector< unsigned > total_vec = list_of< unsigned >
	(1)
	(2)
	(3)
	(10)
	(30)
	;

void
check_calc_ln_gamma_factor_total(const unsigned total)
{
	cout << "******* check_calc_ln_gamma_factor(): total = " << total << "\n";

	obs_vec_t obs(4);

	double uniform = 0.0;
	double dirichlet = 0.0;

	for (obs[0] = 0; obs[0] != total + 1; ++obs[0])
	{
		for (obs[1] = 0; obs[0] + obs[1] != total + 1; ++obs[1])
		{
			for (obs[2] = 0; obs[0] + obs[1] + obs[2] != total + 1; ++obs[2])
			{
				obs[3] = total - obs[0] - obs[1] - obs[2];

				BOOST_ASSERT(unsigned(std::accumulate(obs.begin(), obs.end(), 0)) == total);

				const double ln_uniform = calc_multinomial_ln_likelihood_uniform_dist(obs.begin(), obs.end());
				double p = 0.25;
				BOOST_CHECK_CLOSE(
					ln_uniform,
					calc_multinomial_ln_likelihood(obs.begin(), obs.end(), single_value_iterator< double >( p )),
					0.001);

				uniform += BIO_GSL_EXP(ln_uniform);
				dirichlet += BIO_GSL_EXP(calc_multinomial_ln_likelihood_dirichlet_prior(obs.begin(), obs.end(), alpha.begin()));
			}
		}
	}

	BOOST_CHECK_CLOSE(uniform, 1.0, 0.001);
	BOOST_CHECK_CLOSE(dirichlet, 1.0, 0.001);
}

void
check_calc_ln_gamma_factor(const obs_vec_t & obs)
{
	cout << "******* check_calc_ln_gamma_factor()" << endl;

	boost::io::ios_base_all_saver ias(cout);



	const double n_choose = gsl_sf_exp(calc_ln_n_choose(obs.begin(), obs.end()));
	const double uniform_ll = calc_multinomial_ln_likelihood_uniform_dist(obs.begin(), obs.end());
	const double dirichlet_ll = calc_multinomial_ln_likelihood_dirichlet_prior(obs.begin(), obs.end(), alpha.begin());

	cout.setf(ios::left, ios::adjustfield);
	cout.setf(ios::fixed, ios::floatfield);

#ifdef VERBOSE_CHECKING
	copy(obs.begin(), obs.end(), ostream_iterator< unsigned >(cout , ","));
	cout
		<< " " << setw(6) << uniform_ll
		<< " " << setw(6) << dirichlet_ll
		<< " " << setw(6) << uniform_ll - dirichlet_ll
		<< "\n";
#endif

}



void
register_math_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_PARAM_TEST_CASE( &check_calc_ln_gamma_factor, obs_vec.begin(), obs_vec.end() ), 0);
	test->add( BOOST_PARAM_TEST_CASE( &check_calc_ln_gamma_factor_total, total_vec.begin(), total_vec.end() ), 0);
}
