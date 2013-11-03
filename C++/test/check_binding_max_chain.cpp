/**
@file

Copyright John Reid 2006

*/

#include "bio/defs.h"
#include "bio/bio_max_chain.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;
using namespace std;
USING_BIO_NS;


//BOOST_TEST_DONT_PRINT_LOG_VALUE( max_chain_test::lcs::string )

//#define VERBOSE_CHECKING

namespace check_binding_max_chain {

typedef std::vector< match_result_vec_t > match_results_array;
typedef boost::tuple<
	std::string, 
	match_results_array
> test_case;
typedef std::vector< test_case > test_case_vec;

extern test_case_vec test_cases;


void
check_binding_hit_max_chain( const test_case & tc )
{
	cout << "******* check_binding_hit_max_chain():" << tc.get< 0 >() << "\n";

	binding_hit_max_chain( tc.get< 1 >() );
}


test_case_vec test_cases = tuple_list_of< std::string, match_results_array >
	(
		std::string( "test 2" ),
		list_of< match_result_vec_t >
			(
				list_of< MatchResults >
						( MatchResults( TableLink( MATRIX_DATA, 00001 ), Hit( 0.4f,   3 ) ) )
			)
			(
				list_of< MatchResults >
					( MatchResults( TableLink( MATRIX_DATA, 00001 ), Hit( 0.4f,   4 ) ) )
			)
	)
	(
		std::string( "test 1" ),
		match_results_array()
	)
	;



} //namespace check_binding_max_chain

void
register_binding_max_chain_tests( boost::unit_test::test_suite * test )
{
	using namespace check_binding_max_chain;

	test->add( BOOST_PARAM_TEST_CASE( &check_binding_hit_max_chain, test_cases.begin(), test_cases.end() ) );
}

