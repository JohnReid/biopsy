/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#include "bio_test_defs.h"
#include "bio/equivalent_factors.h"
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



//#define VERBOSE_CHECKING




void
check_factor_synonyms()
{
	cout << "******* check_factor_synonyms()" << endl;

	EquivalentFactors::ptr_t equiv_factors = EquivalentFactors::construct_from_biobase();

	equiv_factors->get_partition( 84 );


#ifdef VERBOSE_CHECKING
#endif

}



void
register_factor_synonym_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_factor_synonyms ), 0);
}


//
// Are we going to compile this test into its own executable?
//
#ifdef BIO_STANDALONE_TEST
BIO_DEFINE_STANDALONE_TEST( "factor synonyms", register_factor_synonym_tests )
#endif //BIO_STANDALONE_TEST
