

#include "bio/score_map.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;
using namespace boost::multi_index;
using boost::multi_index_container;

#include <iostream>
#include <string>
#include <iostream>
#include <iterator>
using namespace std;





//#define VERBOSE_CHECKING


typedef ScoreMap< std::string > score_map_t;



BIO_NS_START //annoyingly need this for std::ostream_iterator to spot this definition...
std::ostream &
operator<<( std::ostream & os, const score_map_t::value_type & value )
{
	os << "(" << value.key << "," << value.score << ")";
	return os;
}
BIO_NS_END




void
check_score_map()
{
	std::cout << "******* check_score_map()\n";

	const score_map_t::type test_map = boost::assign::list_of
		( score_map_t::value_type( "c", 0.3 ) )
		( score_map_t::value_type( "d", 0.4 ) )
		( score_map_t::value_type( "e", 0.5 ) )
		( score_map_t::value_type( "z", 0.4 ) )
		;

#ifdef VERBOSE_CHECKING
	std::copy( test_map.get<0>().begin(), test_map.get<0>().end(), std::ostream_iterator< score_map_t::value_type >( std::cout, "\n" ) );
	std::cout << "\n";

	std::copy( test_map.get<1>().begin(), test_map.get<1>().end(), std::ostream_iterator< score_map_t::value_type >( std::cout, "\n" ) );
	std::cout << "\n";
#endif //VERBOSE_CHECKING

	score_map_t::key_index::const_iterator k = test_map.get<0>().begin();
	BOOST_CHECK_CLOSE((k++)->score, 0.3, 1.0);
	BOOST_CHECK_CLOSE((k++)->score, 0.4, 1.0);
	BOOST_CHECK_CLOSE((k++)->score, 0.5, 1.0);
	BOOST_CHECK_CLOSE((k++)->score, 0.4, 1.0);

	score_map_t::score_index::const_iterator s = test_map.get<1>().begin();
	BOOST_CHECK_CLOSE((s++)->score, 0.3, 1.0);
	BOOST_CHECK_CLOSE((s++)->score, 0.4, 1.0);
	BOOST_CHECK_CLOSE((s++)->score, 0.4, 1.0);
	BOOST_CHECK_CLOSE((s++)->score, 0.5, 1.0);

}

void
register_score_map_tests(boost::unit_test::test_suite * test_map)
{
	test_map->add( BOOST_TEST_CASE( &check_score_map ), 0);
}
