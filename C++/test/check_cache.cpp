

#include "bio/cache.h"
#include "bio/singleton.h"
#include "bio/unary_compose.h"
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



struct CacheElement
{
	typedef boost::shared_ptr< CacheElement > ptr_t;

	int i;

	CacheElement( int i )
		: i( i )
	{
	}

	CacheElement( const CacheElement & rhs )
		: i( rhs.i )
	{
	}

	~CacheElement()
	{
	}
};

struct CacheElementFactory
	: std::unary_function< int, CacheElement::ptr_t >
{
	CacheElement::ptr_t operator()( int i ) const
	{
		return CacheElement::ptr_t( new CacheElement( i ) );
	}
};


struct CacheTest
	: unary_compose<
		Dereference< CacheElement >,
		Cache< CacheElementFactory >
	>
	, Singleton< CacheTest >
{
};



void
check_cache()
{
	cout << "******* check_cache()" << endl;

	CacheTest cache;

	cache(0);
	cache(1);

#ifdef VERBOSE_CHECKING
#endif

}



void
register_cache_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_cache ), 0);
}
