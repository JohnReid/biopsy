
#include <bio/matrix_match.h>
USING_BIO_NS;

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;


#include <iostream>
using namespace std;

void
check_matrix_match_map()
{
	cout << "******* check_matrix_match_map()" << endl;

	get_min_fp_match_map().find( TableLink( MATRIX_DATA, 1 ) );
}


void register_matrix_match_map_tests(test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_matrix_match_map ), 0);
}


