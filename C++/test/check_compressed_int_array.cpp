

#include "bio/compressed_int_array.h"

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
check_compressed_int_array()
{
	cout << "******* check_compressed_int_array()" << endl;

	typedef short value_int;

	const std::vector< value_int > values = list_of
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		(  6 )
		(  0 )
		(  1 )
		(  2 )
		(  3 )
		(  4 )
		(  5 )
		(  6 )
		(  7 )
		;

	typedef compressed_int_array< unsigned, 4 > array;

	BOOST_CHECK_EQUAL( sizeof( unsigned ) * 8, 32 );

	BOOST_CHECK_EQUAL( array::storage_idx(  0 ), 0 );
	BOOST_CHECK_EQUAL( array::storage_idx(  7 ), 0 );
	BOOST_CHECK_EQUAL( array::storage_idx(  8 ), 1 );
	BOOST_CHECK_EQUAL( array::storage_idx( 15 ), 1 );
	BOOST_CHECK_EQUAL( array::storage_idx( 16 ), 2 );

	BOOST_CHECK_EQUAL( array::in_storage_idx(  0 ), 0 );
	BOOST_CHECK_EQUAL( array::in_storage_idx(  7 ), 7 );
	BOOST_CHECK_EQUAL( array::in_storage_idx(  8 ), 0 );
	BOOST_CHECK_EQUAL( array::in_storage_idx( 15 ), 7 );
	BOOST_CHECK_EQUAL( array::in_storage_idx( 16 ), 0 );

	array a( values.size() );
	for( unsigned i = 0; values.size() != i; ++i )
	{
		a.set( i, values[i] );
		BOOST_CHECK_EQUAL( a.get< value_int >( i ), values[i] );
	}
}



void
register_compressed_int_array_tests( boost::unit_test::test_suite * test )
{
	test->add( BOOST_TEST_CASE( &check_compressed_int_array ), 0);
}

test_suite*
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite* suite = BOOST_TEST_SUITE( "Bio test suite" );

	try
	{
		register_compressed_int_array_tests( suite );
	}
	catch (const exception & e)
	{
		cerr << "Exception: " << e.what() << endl;
	}

    return suite;
}

