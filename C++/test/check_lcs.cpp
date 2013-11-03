/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#include "bio_test_defs.h"
#include "bio/defs.h"
#include "bio/lcs.h"
#include "bio/max_chain.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/range/algorithm/equal.hpp>

using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;
using namespace std;
USING_BIO_NS;



//#define VERBOSE_CHECKING


struct random_gen
{
	typedef boost::mt19937 rng;                 // produces randomness out of thin air
	typedef boost::uniform_int< int > int_dist;
	typedef boost::uniform_int< unsigned > unsigned_dist;
	typedef boost::uniform_real< double > weight_dist;
	typedef boost::variate_generator< rng &, int_dist > int_gen;
	typedef boost::variate_generator< rng &, unsigned_dist > unsigned_gen;
	typedef boost::variate_generator< rng &, weight_dist > weight_gen;


	unsigned_dist _num_points_dist;
	unsigned_gen _num_points_gen;

	unsigned_dist _char_dist;
	unsigned_gen _char_gen;

	int_dist _coord_dist;
	int_gen _coord_gen;

	weight_dist _weight_dist;
	weight_gen _weight_gen;

	random_gen( 
		rng & _rng,
		unsigned max_num_points,
		int max_coord,
		unsigned alphabet_size = 26 ) 
		: _num_points_dist( 0, max_num_points )
		, _num_points_gen( _rng, _num_points_dist )
		, _char_dist( 0, alphabet_size - 1 )
		, _char_gen( _rng, _char_dist )
		, _coord_dist( 0, max_coord )
		, _coord_gen( _rng, _coord_dist )
		, _weight_dist( 0.0, 1.0 )
		, _weight_gen( _rng, _weight_dist )
	{
	}

	char character( )
	{
		return 'a' + _char_gen();
	}

	int coord( )
	{
		return _coord_gen();
	}

	std::pair< int, int > coord_range( )
	{
		std::pair< int, int > result( coord(), coord() );
		if( result.first > result.second )
		{
			std::swap( result.first, result.second );
		}
		return result;
	}

	int num_points( )
	{
		return _num_points_gen();
	}

	double weight( )
	{
		return _weight_gen();
	}

	template< unsigned d >
	boost::array< int, d > point()
	{
		typedef boost::array< int, d > point;
		point p;
		std::generate( p.begin(), p.end(), boost::bind< int >( &random_gen::coord, this ) );
		return p;
	}

	template< unsigned d >
	std::pair< boost::array< int, d >, boost::array< int, d > > range()
	{
		typedef std::pair< boost::array< int, d >, boost::array< int, d > > point_pair;
		point_pair points( point< d >(), point< d >() );
		//make sure first <= second in every dimension
		for( unsigned i = 0; d != i; ++i )
		{
			if( points.first[ i ] > points.second[ i ] )
			{
				std::swap( points.first[ i ], points.second[ i ] );
			}
		}
		return points;
	}

};


namespace range_tree {

	
template< unsigned d >
struct test
{
	struct traits : point_traits< int, d >
	{
		typedef double weight;
	};

	typedef std::pair< typename traits::coord_array, typename traits::weight > test_point;
	typedef std::vector< test_point > test_point_vec;

	typedef RangeTree< d, traits > range_tree;
	typedef typename range_tree_ns::node< d, traits, d > node;
	typedef typename node::ptr tree;

	random_gen::rng _rng; //random number generator

	test() : _rng( 42u ) { } //seed the rng

	void generate_and_check( unsigned num_test_cases, unsigned max_points, unsigned max_coord )
	{
		cout 
			<< "******* lcs::range_tree::test::check(): "
			<< num_test_cases << " test cases of up to " 
			<< max_points << " points of dimension " << d << "\n";

		for( unsigned i = 0; num_test_cases != i; ++i )
		{
			check( create_case( max_points, max_coord ) );
		}
	}

	test_point_vec create_case( unsigned max_points, unsigned max_coord )
	{
		random_gen gen( _rng, max_points, max_coord );

		const unsigned num_points = gen.num_points();				// generate a random number of points

		test_point_vec test_case;
		for( unsigned i = 0; num_points != i; ++i )
		{
			test_point p( gen.point< d >(), gen.weight() );
			test_case.push_back( p );
		}

		return test_case;
	}


	void print_point( typename traits::point p )
	{
		for( unsigned i = 0; d != i; ++i )
		{
			std::cout << ( *p )[ i ];
			if( d - 1 != i )
			{
				std::cout << ",";
			}
		}
	}

	void
	check( const test_point_vec & test_points )
	{
		using namespace range_tree;
		using namespace boost::lambda;
		using boost::lambda::_1;

		//get the set of all points...
		typedef std::set< typename traits::point, typename traits::point_less > point_set;
		point_set points;
		for( unsigned i = 0; test_points.size() != i; ++i )
		{
			typename traits::point p = boost::addressof( test_points[ i ].first );
			points.find( p );
			points.insert( p );
#ifdef VERBOSE_CHECKING
			print_point( p );
			std::cout << std::endl;
#endif //VERBOSE_CHECKING
		}
#ifdef VERBOSE_CHECKING
		std::cout << std::endl;
#endif //VERBOSE_CHECKING

		//build tree
		tree t;
		typename range_tree_ns::create_tree< d, traits, d >()( t, points );
		typename range_tree_ns::check_tree< d, traits, d >()( t );
#ifdef VERBOSE_CHECKING
		typename range_tree_ns::visit_tree< d, traits, d >()( t, typename range_tree::print_point< d, traits, d >() );
#endif //VERBOSE_CHECKING

		//check tree has the correct points
		point_set tree_points;
		typename range_tree_ns::visit_tree< d, traits, d >()(
			t,
			range_tree_ns::make_get_points< d, traits, d >(
				std::inserter( tree_points, tree_points.begin() ) ) );
		BOOST_ASSERT( points == tree_points );

		//insert points and weights
		BOOST_FOREACH( const typename test_point_vec::value_type & v, test_points )
		{
			typename range_tree_ns::insert_point< d, traits, d >()( t, boost::addressof( v.first ), v.second );
			//std::cout << v.first[ 0 ] << "," << v.first[ 1 ] << " - " << max_chain::max_weight< 2 >()( t, addressof( v.first ) ) << "\n";
		}
		//std::cout << "\n";

#ifdef VERBOSE_CHECKING
		//query points
		BOOST_FOREACH( traits::point p, points )
		{
			std::copy( p->begin(), p->end(), std::ostream_iterator< traits::coord >( std::cout, "," ) );
			std::cout << " - " << range_tree::max_weight< d >()( t, p ) << "\n";
		}
#endif //VERBOSE_CHECKING

		const typename range_tree::heaviest_point null_point =
			typename range_tree_ns::max_weight< d, traits, d >()( t, traits::always_dominated_point() );
		BOOST_CHECK_EQUAL( null_point.first, traits::always_dominated_point() );
		BOOST_CHECK_EQUAL( null_point.second, 0.0 );

		BOOST_FOREACH( typename traits::point p1, points )
		{
			const typename range_tree::heaviest_point h1 =
				typename range_tree_ns::max_weight< d, traits, d >()( t, p1 );

			BOOST_FOREACH( typename traits::point p2, points )
			{
				if( typename traits::point_dominate()( p1, p2 ) )
				{
					const typename range_tree::heaviest_point h2 =
						range_tree_ns::max_weight< d, traits, d >()( t, p2 );

					BOOST_CHECK( h1.second >= h2.second );
				}
			}
		}
	}
};


void 
generate_and_check( )
{
	const unsigned num_tests = 10;
	const std::vector< unsigned > max_points = list_of( 4 )( 15 )( 300 );

	BOOST_FOREACH( unsigned mp, max_points )
	{
		const unsigned max_coord = mp;

		test<  1 >().generate_and_check( num_tests, mp, max_coord );
		test<  2 >().generate_and_check( num_tests, mp, max_coord );
#ifndef _MSC_VER
		test<  3 >().generate_and_check( num_tests, mp, max_coord );
		test<  4 >().generate_and_check( num_tests, mp, max_coord );

#if 0
		test<  5 >().generate_and_check( num_tests, mp, max_coord );
		test<  6 >().generate_and_check( num_tests, mp, max_coord );
		test<  7 >().generate_and_check( num_tests, mp, max_coord );
		test<  7 >().generate_and_check( num_tests, mp, max_coord );
		test<  8 >().generate_and_check( num_tests, mp, max_coord );
		test<  9 >().generate_and_check( num_tests, mp, max_coord );
#endif
#endif //_MSC_VER
	}
}



}


namespace max_chain_test {

struct V
{
	typedef char character;

	V( char c = ' ', int start = -1, int end = 0, double score = 1.0 )
		: _start( start )
		, _end( end )
		, _c( c )
		, _score( score )
	{
		if( _start >= end )
		{
			throw std::invalid_argument( "start must be < end" );
		}
	}

	int _start;
	int _end;
	char _c;
	double _score;
};

struct VTraits
{
    typedef V hit;
    typedef const V * data;
	typedef char character;
	typedef double weight;
	typedef int coord;

	static data get_data( const V & v ) { return boost::addressof( v ); }
	static character get_char( const V & v ) { return v._c; }
	static weight get_weight( const V & v ) { return v._score; }
	static coord get_start( const V & v ) { return v._start; }
	static coord get_end( const V & v ) { return v._end; }
};


struct CharExtractor
{
	char operator()( const V & v ) const
	{
		return v._c;
	}
};

struct StartExtractor
{
	int operator()( const V & v ) const
	{
		return v._start;
	}
};

struct EndExtractor
{
	int operator()( const V & v ) const
	{
		return v._end;
	}
};

struct ScoreExtractor
{
	double operator()( const V & v ) const
	{
		return v._score;
	}
};

template< unsigned d >
struct VBox
{
	struct traits : point_traits< int, d >
	{
		typedef VBox< d > box;
		typedef const box * box_ptr;
		typedef double weight;
		typedef const V * data;
		typedef std::vector< box > box_vec;					/**< Vector of boxes. */
		typedef std::list< box_ptr > box_ptr_list;			/**< List of box pointers. */
		typedef typename point_traits< int, d >::coord_array coord_array;
		typedef typename point_traits< int, d >::point point;

		static unsigned get_dimensions() { return d; }

		static box make_box( data _d, weight w, coord_array s, coord_array e ) { return box( _d, w, s, e ); }
		static point get_start( const box & b ) { return boost::addressof( b._start ); }
		static point get_end( const box & b ) { return boost::addressof( b._end ); }
	};

	typedef typename traits::box_ptr box_ptr;
	typedef typename traits::weight weight;
	typedef typename traits::data data;
	typedef typename traits::coord_array point;

	VBox( 
		typename traits::data _d, 
		typename traits::weight w, 
		typename traits::coord_array s, 
		typename traits::coord_array e )
		: _data( _d )
		, _weight( w )
		, _start( s )
		, _end( e )
	{
	}

    typename traits::data _data;
	typename traits::weight _weight;
	typename traits::coord_array _start;
	typename traits::coord_array _end;

};




typedef LCS< V, char, CharExtractor, StartExtractor, EndExtractor, ScoreExtractor > lcs;

typedef std::vector< V > v_vec;
typedef std::vector< v_vec > v_array;
typedef boost::tuple< std::string, v_array, lcs::string, double > lcs_test_case;
typedef std::vector< lcs_test_case > lcs_test_case_vec;

} // namespace max_chain_test


bool
operator==( const max_chain_test::lcs::string & lhs, const max_chain_test::lcs::string & rhs )
{
    return boost::equal( lhs, rhs );
}

std::ostream &
operator<<( std::ostream & os, const max_chain_test::lcs::string & s )
{
    os << "\"";
    BOOST_FOREACH( max_chain_test::lcs::character c, s )
    {
        os << c;
    }
    os << "\"";
    return os;
}




namespace max_chain_test {

struct generate_lcs_test_case
{
	random_gen::rng _rng;
	generate_lcs_test_case() : _rng( 42u ) { }

	template< typename Int >
	Int gen_int( Int max )
	{
		typedef uniform_int< Int > dist;
		return boost::variate_generator< random_gen::rng &, dist >( _rng, dist( 0, max ) )();
	};

	template< typename Real >
	Real gen_real( )
	{
		typedef uniform_real< Real > dist;
		return boost::variate_generator< random_gen::rng &, dist >( _rng, dist( 0.0, 1.0 ) )();
	};

	v_array operator()( unsigned num_sequences, unsigned max_sequence_length, unsigned max_coord, unsigned alphabet_size )
	{
		v_array result( num_sequences );
		for( unsigned i = 0; num_sequences != i; ++i )
		{
			const unsigned length = gen_int< unsigned >( max_sequence_length );
			result[ i ].resize( length );
			for( unsigned j = 0; length != j; ++j )
			{
				const char c = 'a' + gen_int< unsigned >( alphabet_size - 1 );
				int start = gen_int< int >( max_coord );
				int end = gen_int< int >( max_coord );
				if( start > end )
				{
					std::swap( start, end );
				}
				else if( start == end )
				{
					++end;
				}
				const double score = gen_real< double >();
				result[ i ][ j ] = V( c, start, end, score );
			}
		}
		return result;
	}
};


//using boost::assign::tuple_list_of;
const lcs_test_case_vec test_cases = tuple_list_of
    (
        "empty"
        ,
        v_array()
        ,
        lcs::string()
        ,
        0.0
    )
    (
        "singular - 2D"
        ,
        list_of< v_vec >
            ( list_of ( V( 'a',  0, 10, 1. ) ) )
            ( list_of ( V( 'a', 10, 20, 2. ) ) )
        ,
        list_of( 'a' )
        ,
        1.5
    )
    (
        "singular - 3D"
        ,
        list_of< v_vec >
            ( list_of ( V( 'a',  0, 10, 1. ) ) )
            ( list_of ( V( 'a', 10, 20, 2. ) ) )
            ( list_of ( V( 'a', 20, 30, 3. ) ) )
        ,
        list_of( 'a' )
        ,
        2.
    )
    (
        "singular - 4D"
        ,
        list_of< v_vec >
            ( list_of ( V( 'a',  0, 10, 1. ) ) )
            ( list_of ( V( 'a', 10, 20, 2. ) ) )
            ( list_of ( V( 'a', 20, 30, 3. ) ) )
            ( list_of ( V( 'a', 30, 40, 4. ) ) )
        ,
        list_of( 'a' )
        ,
        2.5
    )
	(
		"problem from random tests - max_chain swaps first two matches "
		,
		list_of< v_vec >
			( list_of
                ( V( 'd', 17, 20, 0.996334 ) )( V( 'f', 6, 14,  0.448611 ) )( V( 'h', 9, 14, 0.886196 ) )
                ( V( 'e', 7, 16,  0.391526 ) )( V( 'h', 15, 20, 0.695618 ) )( V( 'f', 5, 7, 0.619589 ) )
                ( V( 'b', 14, 15, 0.794197 ) )( V( 'd', 5, 19,  0.588202 ) )( V( 'a', 13, 20, 0.642326 ) )
                ( V( 'b', 6, 13,  0.579984 ) )( V( 'b', 0, 4,    0.56066 ) )( V( 'i', 8, 17, 0.676468 ) )
                ( V( 'h', 0, 18,  0.269821 ) )( V( 'g', 13, 17, 0.498256 ) )( V( 'd', 1, 2, 0.0585509 ) )
                ( V( 'e', 8, 20,  0.784897 ) ) )
			( list_of
                ( V( 'j', 13, 18, 0.825131 ) )( V( 'j', 15, 20, 0.373463 ) )( V( 'b', 0, 15, 0.735506 ) )
                ( V( 'i', 0, 11,  0.371422 ) )( V( 'b', 7, 8,   0.801367 ) )( V('g', 0, 9, 0.733616 ) )
                ( V( 'b', 13, 17,  0.17348 ) )( V( 'c', 3, 14,  0.496166 ) )( V( 'j', 3, 11, 0.685911 ) )
                ( V( 'j', 12, 13,  0.74722 ) )( V( 'g', 1, 13,  0.950313 ) )( V( 'd', 3, 20, 0.56555 ) )
                ( V( 'b', 3, 10,  0.165294 ) )( V( 'd', 3, 9,    0.58581 ) )( V( 'd', 10, 16, 0.697671 ) )
                ( V( 'b', 11, 13, 0.782582 ) )( V( 'g', 3, 18,  0.852728 ) )( V( 'j', 9, 20, 0.376914 ) )
                ( V( 'e', 2, 12,  0.480906 ) )( V( 'd', 15, 16, 0.937353 ) )( V( 'j', 1, 8, 0.631171 ) )
                ( V( 'f', 9, 20,  0.903913 ) )( V( 'g', 8, 15,  0.780213 ) )( V( 'e', 13, 15, 0.0572803 ) )
                ( V( 'f', 0, 18,  0.131349 ) )( V( 'h', 2, 12,  0.922623 ) )( V( 'f', 3, 4, 0.0134026 ) )
                ( V( 'b', 5, 6,  0.0445949 ) )( V( 'h', 1, 3,    0.27572 ) )( V( 'f', 10, 11, 0.787127 ) ) )
		,
		list_of( 'b' )( 'f' )( 'b' )( 'd' )
		,
		3.1396045
	)
	(
		"problem with start == end in first V"
		,
		list_of< v_vec >
			( list_of( V( 'b', 5, 6, 0.386 ) )( V( 'f', 3, 11, 0.191 ) )( V( 'h', 4, 11, 0.619 ) )( V( 'b', 12, 14, 0.484 ) )( V( 'c', 14, 18, 0.602 ) ) )
		,
		list_of( 'h' )( 'b' )( 'c' )
		,
		1.705
	)
	(
		"one good hit score"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 1, 0.1 ) )( V( 'b', 1, 2, 0.1 ) )( V( 'c', 2, 3, 0.1 ) )( V( 'd', 3, 4, 0.1 ) ) )
			( list_of( V( 'd', 0, 1, 0.9 ) )( V( 'a', 0, 1, 0.1 ) )( V( 'b', 1, 2, 0.1 ) )( V( 'c', 2, 3, 0.1 ) )( V( 'd', 3, 4, 0.1 ) ) )
		,
		list_of( 'd' )
		,
		0.5
	)
	(
		"long chain score"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 1, 0.1 ) )( V( 'b', 1, 2, 0.1 ) )( V( 'c', 2, 3, 0.1 ) )( V( 'd', 3, 4, 0.1 ) ) )
			( list_of( V( 'd', 0, 1, 0.3 ) )( V( 'a', 0, 1, 0.1 ) )( V( 'b', 1, 2, 0.1 ) )( V( 'c', 2, 3, 0.1 ) )( V( 'd', 3, 4, 0.1 ) ) )
		,
		list_of( 'a' )( 'b' )( 'c' )( 'd' )
		,
		0.4
	)
	(
		"next to"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
			( list_of( V( 'a', 0, 5 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) ) )
		,
		list_of( 'c' )( 'g' )
		,
		2.0
	)
	(
		"overlapping"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 3 ) )( V( 'c', 1, 3 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
			( list_of( V( 'a', 0, 9 ) )( V( 'c', 1, 9 ) )( V( 'g', 2, 9 ) ) )
		,
		list_of( 'a' )
		,
		1.0
	)
	(
		"empty sequence"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
			( v_vec() )
		,
		lcs::string()
		,
		0.0
	)
	(
		"empty sequence 2"
		,
		list_of< v_vec >
			( v_vec() )
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
		,
		lcs::string()
		,
		0.0
	)
	(
		"identical"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
			( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
		,
		list_of( 'a' )( 'c' )( 'g' )( 't' )
		,
		4.0
	)
	(
		"score"
		,
		list_of< v_vec >
			( list_of( V( 'a', 0, 1, 0.1 ) )( V( 'c', 1, 2, 1.0 ) ) )
			( list_of( V( 'c', 0, 1, 1.0 ) )( V( 'c', 1, 2, 0.1 ) ) )
		,
		list_of( 'c' )
		,
		1.0
	)
	;


struct int_cmp
{
	bool operator()( int lhs, int rhs ) { return lhs < rhs; }
};

typedef std::pair< lcs::string, double > lcs_result;

template< unsigned d >
lcs_result
get_lcs_using_max_chain_d( const v_array & sequences, bool print = false )
{
    typedef typename VBox< d >::traits box_traits;
    typename box_traits::box_vec boxes;
    auto values = box_generator< VTraits >::template boxes_from_sequences< box_traits >(
        sequences,
        std::back_inserter( boxes ),
        50000
    );

	typedef typename VBox< d >::traits traits;
	typedef MaxChain< d, VBox< d >, traits > max_chain;
	max_chain _max_chain( boxes );

	lcs::string _longest;
	double weight = 0.0;
	BOOST_FOREACH( typename max_chain::box_ptr b, _max_chain.maximal_chain )
	{
		_longest.push_back( b->_data->_c );
#if 0
		std::cout 
			<< b->_data->_start << ","
			<< b->_data->_end << ","
			<< b->_data->_c << ","
			<< b->_data->_score << endl;
#endif
		weight += b->_weight;
	}
	BOOST_CHECK_CLOSE( _max_chain.maximal_weight, weight, 0.001 );

	return std::make_pair( _longest, _max_chain.maximal_weight );
}

lcs_result
get_lcs_using_max_chain( const v_array & sequences, bool print = false )
{
	switch( boost::size( sequences ) )
	{
	case 0: return get_lcs_using_max_chain_d< 0 >( sequences, print ); break;
	case 1: return get_lcs_using_max_chain_d< 1 >( sequences, print ); break;
	case 2: return get_lcs_using_max_chain_d< 2 >( sequences, print ); break;
	case 3: return get_lcs_using_max_chain_d< 3 >( sequences, print ); break;
	case 4: return get_lcs_using_max_chain_d< 4 >( sequences, print ); break;
	case 5: return get_lcs_using_max_chain_d< 5 >( sequences, print ); break;
	case 6: return get_lcs_using_max_chain_d< 6 >( sequences, print ); break;
	case 7: return get_lcs_using_max_chain_d< 7 >( sequences, print ); break;
	//case 8: return get_lcs_using_max_chain_d< 8 >( sequences, print ); break;
	default: throw std::logic_error( "Too many test input sequences" );
	}
}

void
check_max_chain( const lcs_test_case & test_case )
{
    using namespace boost;

	cout << "******* check_max_chain(): \"" << test_case.get< 0 >() << "\"\n";

	if( 1 != test_case.get< 1 >().size() ) //does not work with one sequence at the moment
	{
		lcs_result result = get_lcs_using_max_chain( test_case.get< 1 >() );
//        std::cout
//            << size( test_case.get< 2 >() ) << " : " << size( result.first ) << ", "
//            << test_case.get< 2 >() << " : " << result.first
//            << ", weights " << test_case.get< 3 >() << " : " << result.second << "\n";

//        if( ! boost::equal( test_case.get< 2 >(), result.first ) ) {
//            std::cout << "Not equal under boost range equality\n";
//        }
//        if( test_case.get< 2 >() != result.first) {
//            std::cout << "Not equal under != operator\n";
//        }
		BOOST_CHECK_EQUAL( test_case.get< 2 >(), result.first );
		BOOST_CHECK_CLOSE( test_case.get< 3 >(), result.second, 0.001 );
	}
}

lcs_result
get_lcs_using_lcs_algorithm( const v_array & sequences )
{
	lcs _lcs( sequences );
	_lcs.calculate_best();

	return lcs_result( _lcs.get_best().get_string(), _lcs.get_best().get_score() );
}

void
check_random_lcs_test_cases( )
{
	typedef boost::array< unsigned, 5 > test_params;
	typedef std::vector< test_params > test_params_vec;

	// # tests, # seqs, max seq length, alphabet_size, max coord
	const test_params_vec _params = list_of
		( list_of( 300 )( 0 )( 100 )( 20 )( 100 ) )
		( list_of(   1 )( 1 )( 100 )( 20 )( 100 ) )
		( list_of( 300 )( 2 )( 100 )( 20 )( 100 ) )
		( list_of( 300 )( 3 )( 100 )( 20 )( 100 ) )
		( list_of( 300 )( 4 )(  30 )( 20 )( 100 ) )
		( list_of( 300 )( 5 )(  30 )( 20 )( 100 ) )
		( list_of( 300 )( 5 )(  30 )(  4 )( 100 ) )
		( list_of( 300 )( 5 )(  30 )(  4 )(  10 ) )
		( list_of( 300 )( 6 )(  30 )(  4 )(  10 ) )
		( list_of( 300 )( 7 )(  30 )(  4 )(  10 ) )
		;

	generate_lcs_test_case gen;

	BOOST_FOREACH( const test_params & p, _params )
	{
		cout 
			<< "******* check_random_lcs_test_cases(): checking " << p[ 0 ] << " random cases\n"
			<< "# seqs: " << p[ 1 ] << "\n"
			<< "max length: " << p[ 2 ] << "\n"
			<< "alphabet size: " << p[ 3 ] << "\n"
			<< "max coord: " << p[ 4 ] << "\n"
			;

		double lcs_cumulative_time = 0.0;
		double mc_cumulative_time = 0.0;
		for( unsigned i = 0; p[ 0 ] != i; ++i )
		{
			//generate a test case
			v_array sequences( gen( p[ 1 ], p[ 2 ], p[ 3 ], p[ 4 ] ) );
#ifdef VERBOSE_CHECKING
			cout << "# seqs: " << sequences.size();
#endif //VERBOSE_CHECKING

			boost::timer _timer;
			const lcs_result _lcs_result = get_lcs_using_lcs_algorithm( sequences );
			lcs_cumulative_time += _timer.elapsed();
#ifdef VERBOSE_CHECKING
			cout  << ", LCS: " << _timer.elapsed();
#endif //VERBOSE_CHECKING
			
			if( 1 == sequences.size() )
			{
				//max_chain does not work for single sequences...
#ifdef VERBOSE_CHECKING
				cout << "\n";
#endif //VERBOSE_CHECKING
				continue;
			}
			
			_timer.restart();
			const lcs_result max_chain_result = get_lcs_using_max_chain( sequences );
			mc_cumulative_time += _timer.elapsed();
#ifdef VERBOSE_CHECKING
			cout  << ", MC: " << _timer.elapsed() << "\n";
#endif //VERBOSE_CHECKING
			
			//do we have a mismatch?
			if( max_chain_result.first != _lcs_result.first )
			{
				get_lcs_using_max_chain( sequences, true ); //run again with print = true
				BOOST_FOREACH( const v_vec & seq, sequences )
				{
					//( list_of( V( 'a', 0, 1 ) )( V( 'c', 1, 2 ) )( V( 'g', 2, 3 ) )( V( 't', 3, 4 ) ) )
					cout << "( list_of";
					BOOST_FOREACH( const V & v, seq )
					{
						cout << "( V( '" << v._c << "', " << v._start << ", " << v._end << ", " << v._score << " ) )";
					}
					cout << " )\n";
				}
				//throw std::logic_error( "Results don't agree" );
			}

			//check
			BOOST_CHECK_EQUAL( max_chain_result.first, _lcs_result.first );
			BOOST_CHECK_CLOSE( max_chain_result.second, _lcs_result.second, 0.001 );
		}
		std::cout << "LCS:" << lcs_cumulative_time << " secs\n";
		std::cout << "MC :" << mc_cumulative_time << " secs\n";
	}
}

void
check_lcs( const lcs_test_case & test_case )
{
	cout << "******* check_lcs(): \"" << test_case.get< 0 >() << "\"\n";

	lcs _lcs( test_case.get< 1 >() );
	_lcs.calculate_best();

#ifdef VERBOSE_CHECKING
	std::cout << "Universe:\n";
	BOOST_FOREACH( lcs::character c, _lcs._universe )
	{
		std::cout << c << "\n";
	}

	std::cout << "Sequences:\n";
	BOOST_FOREACH( lcs::seq s, _lcs._sequences )
	{
		BOOST_FOREACH( lcs::value c, s )
		{
			std::cout << _lcs.get_char( c  )<< " ";
		}
		std::cout << "\n";
	}

	std::cout << "End positions:\n";
	BOOST_FOREACH( lcs::seq s, _lcs._sequences )
	{
		BOOST_FOREACH( lcs::value v, s )
		{
			std::cout << _lcs.get_end( v  )<< " ";
		}
		std::cout << "\n";
	}

	std::cout << "End position\n";
	BOOST_FOREACH( unsigned i, _lcs._end_position )
	{
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::cout << "End value\n";
	BOOST_FOREACH( unsigned i, _lcs._end_value )
	{
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::cout << "Indices\n";
	lcs::multi_iterator end_it = _lcs.get_end_iterator();
	for( lcs::multi_iterator end_it = _lcs.get_end_iterator(); ! end_it.at_end(); end_it.next() )
	{
		BOOST_FOREACH( unsigned i, end_it._ind )
		{
			std::cout << i << " ";
		}
		std::cout << "\n";

		for( lcs::multi_iterator value_it = _lcs.get_value_range_for_end( end_it._ind ); ! value_it.at_end(); value_it.next() )
		{
			for( unsigned d = 0; _lcs.get_num_seqs() != d; ++d )
			{
				std::cout << " " << _lcs.get_end( _lcs.get_value( d, value_it._ind[ d ] ) );
			}
			std::cout << "\n";
		}
	}
#endif

	BOOST_CHECK_EQUAL( test_case.get< 2 >(), _lcs.get_best().get_string() );
	BOOST_CHECK_CLOSE( test_case.get< 3 >(), _lcs.get_best().get_score(), 0.001 );

	lcs::value_holder::vec values = _lcs.get_best().get_values();

#ifdef VERBOSE_CHECKING
	BOOST_FOREACH( const lcs::value_holder & v, values )
	{
		std::cout
			<< v._start
			<< "," << v._end
			<< "," << v._score
			<< "," << v._c
			<< "\n";
	}
#endif //VERBOSE_CHECKING
}

}

BOOST_TEST_DONT_PRINT_LOG_VALUE( max_chain_test::lcs::string )


#ifdef TEST_AGAINST_CGAL
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

void
check_cgal_range_tree()
{
	cout << "******* check_cgal_range_tree()\n";

	typedef CGAL::Cartesian< int > K;
	typedef CGAL::Range_segment_tree_set_traits_2< K > Traits;
	typedef CGAL::Range_tree_2< Traits > Range_tree_2_type;

	typedef Traits::Key Key;
	typedef std::vector<Key> KeyVec;
	typedef Traits::Interval Interval;
	std::vector<Key>::iterator first, last, current;

	KeyVec InputList = list_of
		( Key( 1, 1) )
		( Key( 1, 2) )
		( Key( 1, 3) )
		( Key( 2, 1) )
		( Key( 2, 2) )
		( Key( 2, 3) )
		( Key( 3, 1) )
		( Key( 3, 2) )
		( Key( 3, 3) )
		;

	//create a range tree
	Range_tree_2_type Range_tree_2( InputList.begin(), InputList.end() );

	const Key l( 1, 1 );
	const Key u( 3, 3 );
	const Interval win = Interval( l, u );
	std::cout
		<< "Window Query: lower left point: "
		<< l
		<< "\nupper right point: "
		<< u
		<< "\n";

	KeyVec OutputList;
	Range_tree_2.window_query( win, std::back_inserter(OutputList) );
	current = OutputList.begin();
	while( current != OutputList.end())
	{
		std::cout << (*current).x()<< "-" << (*current).y() << std::endl;
		current++;
	}
}

#endif //TEST_AGAINST_CGAL






void
register_lcs_tests( boost::unit_test::test_suite * test )
{
	using namespace max_chain_test;
	using namespace range_tree;

    test->add( BOOST_TEST_CASE( &check_random_lcs_test_cases ) );
	test->add( BOOST_PARAM_TEST_CASE( &check_lcs, test_cases.begin(), test_cases.end() ) );
	test->add( BOOST_PARAM_TEST_CASE( &check_max_chain, test_cases.begin(), test_cases.end() ) );
    test->add( BOOST_TEST_CASE( &range_tree::generate_and_check ) );
#ifdef TEST_AGAINST_CGAL
	test->add( BOOST_TEST_CASE( &check_cgal_range_tree ) );
#endif //TEST_AGAINST_CGAL
}



//
// Are we going to compile this test into its own executable?
//
#ifdef BIO_STANDALONE_TEST
BIO_DEFINE_STANDALONE_TEST( "lcs", register_lcs_tests )
#endif //BIO_STANDALONE_TEST
