
#include "bio_test_data.h"

#include <bio/pssm_likelihood.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
#include <bio/pssm_likelihood_cache.h>
#include <bio/environment.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/progress.hpp>
using namespace boost;
using namespace boost::unit_test;

#include <iostream>
#include <numeric>
using namespace std;


//#define VERBOSE_CHECKING


void
check_pssm_likelihood(BiobaseParams const & params)
{
	ensure_biobase_params_built();

	cout << "******* check_pssm_likelihood(): Testing with " << params.name.c_str() << " parameters" << endl;

#ifdef VERBOSE_CHECKING
	progress_timer timer;
	const bool verbose = true;
#else
	const bool verbose = false;
#endif

	BiobaseLikelihoods likelihoods(10);

	//build and calculate the likelihoods of the biobase pssm matching matrix
	pssm_likelihood(
		*params.pssm,
		50000,
		likelihoods);

	//build and calculate the likelihoods of the biobase iupac matching matrix
	pssm_likelihood(
		make_pssm_from_iupac(params.iupac->begin(), params.iupac->end()),
		50000,
		likelihoods);

#ifdef VERBOSE_CHECKING
	std::copy(likelihoods.begin(), likelihoods.end(), ostream_iterator<BIO_NS::float_t>(cout, ","));
	cout << endl;
#endif

	//check the likelihoods add to 1
	BOOST_CHECK_EQUAL(1.0, accumulate(likelihoods.begin(), likelihoods.end(), 0.0));

#ifdef VERBOSE_CHECKING
	const BiobaseLikelihoods * l = PssmLikelihoodCache::singleton().get_likelihoods(params.name);
	std::copy(l->begin(), l->end(), ostream_iterator<BIO_NS::float_t>(cout, ","));
	cout << endl;
#endif

}



struct check_likelihoods
{
	void operator()( BIO_NS::float_t likelihood ) const
	{
		BOOST_CHECK( likelihood >= BIO_NS::float_t( 0.0 ) );
		BOOST_CHECK( likelihood <= BIO_NS::float_t( 1.0 ) );

		if ( likelihood > BIO_NS::float_t( 1.0 ) )
		{
			cout << "Problem\n";
		}
	}
};



void check_pssm_likelihood_2()
{
	cout << "******* check_pssm_likelihood_2()" << endl;

	{
		Matrix * pssm = BiobaseDb::singleton().get_entry< MATRIX_DATA >( 328 );
		BiobaseLikelihoods pssm_likelihoods(BioEnvironment::singleton().num_normalisation_quanta);
		pssm_likelihood(
			make_pssm( pssm ),
			BioEnvironment::singleton().max_pssm_likelihood_map_size,
			pssm_likelihoods);

		std::for_each(
			pssm_likelihoods.begin(),
			pssm_likelihoods.end(),
			check_likelihoods() );
	}

	{
		Site * pssm = BiobaseDb::singleton().get_entry< SITE_DATA >( 4538 );
		BiobaseLikelihoods pssm_likelihoods(BioEnvironment::singleton().num_normalisation_quanta);
		pssm_likelihood(
			make_pssm( pssm ),
			BioEnvironment::singleton().max_pssm_likelihood_map_size,
			pssm_likelihoods);

		std::for_each(
			pssm_likelihoods.begin(),
			pssm_likelihoods.end(),
			check_likelihoods() );
	}

}




void register_pssm_likelihoods_tests(test_suite * test)
{
	ensure_biobase_params_built();

	test->add(BOOST_TEST_CASE(&check_pssm_likelihood_2), 0);
	test->add(BOOST_PARAM_TEST_CASE(&check_pssm_likelihood, biobase_params.begin(), biobase_params.end()), 0);
}


