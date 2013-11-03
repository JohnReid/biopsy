/**
@file

Copyright John Reid 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013

*/

#ifndef BIO_MAX_CHAIN_H_
#define BIO_MAX_CHAIN_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "bio/defs.h"
#include "bio/max_chain_boxes.h"
#include "bio/log.h"

#include <boost/array.hpp>
#include <boost/preprocessor/iteration/local.hpp>


BIO_NS_START

namespace range_tree_ns {



/** A node in a k-dimensional tree. */
template<
	unsigned d,
	typename traits,
	typename node
>
struct base_node
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;

	typedef boost::scoped_ptr< node > ptr;

	ptr _left;											/**< Left child. */
	ptr _right;											/**< Right child. */
	point _point;										/**< Only for leaf nodes. */
	weight _max_weight;									/**< Maximum weight of any sub-node. */
	point _min;											/**< Minimum point of all children. */
	point _max;											/**< Maximum point of all children. */

	base_node( point p = 0, weight w = weight( 0.0 ) )
		: _point( p )
		, _max_weight( w )
		, _min( 0 )
		, _max( 0 )
	{
	}
};

/**
A range tree node in a k-dimensional problem.
A 0-dimensional problem has no data.
A 1-dimensional problem has a tree with no sub-trees.
A 2 or greater-dimensional problem has a tree with sub-trees.
We specialise the 1-dimensional case. The 0-dimensional case is never instantiated.
See http://www.scs.carleton.ca/~michiel/slack.pdf for the algorithm.
*/
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
struct node
	: base_node< d, traits, node< d, traits, k > >
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;

	typedef typename node< d, traits, k >::ptr ptr;
	typedef typename node< d, traits, k - 1 >::ptr sub_tree;

	node( point p = 0 ) : base_node< d, traits, node< d, traits, k > >( p ) { }
	sub_tree _sub_tree;
};

/** Specialise a node in a 1-dimensional tree. We have no sub-tree in this case. */
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct node< d, traits, 1 >
	: base_node< d, traits, node< d, traits, 1 > >
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;

	typedef typename node< d, traits, 1 >::ptr ptr;

	node( point p = 0 ) : base_node< d, traits, node< d, traits, 1 > >( p ) { }
};

/** Apply a functor to the sub-tree. */
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
struct
apply_to_sub_tree
{
	template<
		typename F
	>
	void operator()( typename node< d, traits, k >::ptr & t, F f ) const
	{
		f( t->_sub_tree );
	}
};

/** Apply a functor to the sub-tree. */
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct
apply_to_sub_tree< d, traits, 1 >
{
	template<
		typename F
	>
	void operator()( typename node< d, traits, 1 >::ptr & t, F f ) const
	{
		//nothing to do
	}
};

/** Apply a functor to the sub-tree. Version with return type and default result for no-sub-tree
specialisation. */
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k,
	typename Ret
>
struct
apply_to_sub_tree_ret
{
	Ret _default_result;
	apply_to_sub_tree_ret( Ret default_result ) : _default_result( default_result ) { }

	template<
		typename F
	>
	Ret operator()( typename node< d, traits, k >::ptr & t, F f ) const
	{
		return f( t->_sub_tree );
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	typename Ret
>
struct
apply_to_sub_tree_ret< d, traits, 1, Ret >
{
	Ret _default_result;
	apply_to_sub_tree_ret( Ret default_result ) : _default_result( default_result ) { }

	template<
		typename F
	>
	Ret operator()( typename node< d, traits, 1 >::ptr & t, F f ) const
	{
		return _default_result;
	}
};



template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
bool
is_in_left_child( const typename node< d, traits, k >::ptr & t, typename traits::point p )
{
	BOOST_ASSERT( t->_right );
	BOOST_ASSERT( t->_left );
	return typename traits::point_less_k( k - 1 )( p, t->_right->_min ); //is p is less than smallest node on right?
}


template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
typename node< d, traits, k >::ptr &
get_child_for( const typename node< d, traits, k >::ptr & t, typename traits::point p )
{
	return
		is_in_left_child< d, traits, k >( t, p )
			? t->_left
			: t->_right;
}


template<
	unsigned d,
	typename traits,
	unsigned k
>
struct print_point
{
	void operator()( const typename node< d, traits, k >::ptr & t )
	{
		if( t->_point )
		{
			print( t->_point );
		}
	}
	void print( typename traits::point p )
	{
		std::copy( p->begin(), p->begin() + k, std::ostream_iterator< int >( std::cout, "," ) );
		std::cout << "\n";
	}
};

/** Builds a set of points in the tree. */
template<
	unsigned d,
	typename traits,
	unsigned k,
	typename OutIt
>
struct get_points
{
	OutIt _output_it;
	get_points( OutIt output_it ) : _output_it( output_it ) { }

	void operator()( const typename node< d, traits, k >::ptr & t )
	{
		if( t->_point )
		{
			*_output_it = t->_point;
			++_output_it;
		}
	}
};

template<
	unsigned d,
	typename traits,
	unsigned k,
	typename OutIt
>
get_points< d, traits, k, OutIt >
make_get_points( OutIt output_it )
{
	return get_points< d, traits, k, OutIt >( output_it );
}

template<
	unsigned d,
	typename traits,
	unsigned k
>
struct visit_tree
{
	template< typename Visitor >
	void operator()( typename node< d, traits, k >::ptr & t, Visitor v ) const
	{
		if( t )
		{
			operator()( t->_left, v );
			v( t );
			operator()( t->_right, v );
		}
	}
};

template<
	unsigned d,
	typename traits,
	unsigned k
>
struct sort_points
{
	template<
		typename PointsRange,
		typename PointsVec
	>
	void operator()(
		const PointsRange & points,
		PointsVec & sorted_points ) const
	{
		sorted_points.clear( );

		BOOST_FOREACH( typename traits::point p, points )
		{
			sorted_points.push_back( p );
		}
		std::sort(
			sorted_points.begin(),
			sorted_points.end(),
			typename traits::point_less_k( k - 1 ) );
	}
};

template<
	unsigned d,       /**< Dimension of the range tree. */
	typename traits, /**< The traits of the problem. */
	unsigned k
>
struct max_weight
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef typename traits::point_dominate point_dominate;
	typedef std::pair< point, weight > heaviest_point;

	/**
	Returns point_traits::always_dominated_point() if no point in tree dominated by p.
	*/
	heaviest_point operator()( typename node< d, traits, k >::ptr & t, point p ) const
	{
		heaviest_point result( traits::always_dominated_point(), weight( 0.0 ) );

		//do we have a tree?
		if( ! t )
		{
			//no
		}
		//is it a leaf node?
		else if( t->_point )
		{
			//yes

			//need to check our point dominates it or is the same
			if( point_dominate()( p, t->_point ) )
			{
				//it does
				result.first = t->_point;
				result.second = t->_max_weight;
			}
		}
		else
		{
			//check we have left and right children
			BOOST_ASSERT( t->_left && t->_right );

			//descend the tree looking for our point...
			//should we descend the left or the right child?
			if( is_in_left_child< d, traits, k >( t, p ) )
			{
				//descend left child
				result = operator()( t->_left, p );
			}
			else
			{
				//descend right child - have to check left child for maximal weight
				//If k == 1 there is no sub-tree and we should just check the left
				//child. Similarly if the left child is a leaf, then there is no
				//sub-tree and we may as well just check the node.
				heaviest_point sub_tree_result;
				if( 1 == k || t->_left->_point )
				{
					sub_tree_result = operator()( t->_left, p );
				}
				else
				{
					//check sub-tree
					using namespace boost::lambda;
					using boost::lambda::_1;
					sub_tree_result =
						apply_to_sub_tree_ret< d, traits, k, heaviest_point >( result )(
							t->_left,
							boost::lambda::bind< heaviest_point >(
								max_weight< d, traits, k - 1 >(),
								boost::lambda::_1,
								p ) );
				}
				heaviest_point child_result = operator()( t->_right, p );

				result = child_result.second > sub_tree_result.second ? child_result : sub_tree_result;
			}
		}
		return result;
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct max_weight< d, traits, 0 >
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef std::pair< point, weight > heaviest_point;

	template<
		typename Dummy
	>
	heaviest_point operator()( Dummy & t, point p ) const
	{
		return 0.0;
	}
};


template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
struct insert_point
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef std::pair< point, weight > heaviest_point;

	void operator()( typename node< d, traits, k >::ptr & t, point p, weight w ) const
	{
		//update if new maximum
		if( w > t->_max_weight )
		{
			t->_max_weight = w;
		}

		//is this a leaf node?
		if( t->_point )
		{
			//yes - so the points must be equal
			BOOST_ASSERT( typename traits::point_equal()( t->_point, p ) );
		}
		else
		{
			//check we have left and right children
			BOOST_ASSERT( t->_left && t->_right );

			//should we look in the left or the right child
			typename node< d, traits, k >::ptr & child = get_child_for< d, traits, k >( t, p );
			operator()( child, p, w );

			//also update sub-tree
			using namespace boost::lambda;
			using boost::lambda::_1;
			apply_to_sub_tree< d, traits, k >()(
				t,
				boost::lambda::bind< void >(
					insert_point< d, traits, k - 1 >(),
					boost::lambda::_1,
					p,
					w ) );
		}
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct insert_point< d, traits, 0 >
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef std::pair< point, weight > heaviest_point;

	template<
		typename Dummy
	>
	void operator()( Dummy & t, point p, weight w ) const
	{
	}
};

/**
Creates a tree for k-dimensions of data.
*/
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k					/**< k is the tree dimension. */
>
struct create_tree
{
	typedef typename node< d, traits, k >::ptr tree;

	template<
		typename PointsRange
	>
	void operator()(
		tree & t,
		const PointsRange & points ) const
	{
		std::vector< typename traits::point > sorted_points;
		sort_points< d, traits, k >()( points, sorted_points );
		//we can't have any duplicate elements...
		BOOST_ASSERT( std::adjacent_find( sorted_points.begin(), sorted_points.end(), typename traits::point_equal() ) == sorted_points.end() );

		//build our tree
		operator()(
			t,
			sorted_points,
			0,
			points.size() );
	}

	template<
		typename PointsRange
	>
	void operator()(
		tree & t,
		const PointsRange & sorted_points,
		unsigned begin,
		unsigned end ) const
	{
		if( begin == end )
		{
			return;
		}

		BOOST_ASSERT( 0 == t ); //why would the tree already be populated?
		BOOST_ASSERT( begin < end );

		//have we finished our recursion? I.e. is this a leaf node?
		if( begin + 1 == end )
		{
			typename traits::point p = sorted_points[ begin ];
			t.reset( new typename tree::element_type( p ) );
			t->_max = t->_min = p;
			return;
		}
		else //it is an internal node
		{
			//split begin and end in two
			const unsigned middle = ( begin + end ) / 2;
			BOOST_ASSERT( middle != begin );
			BOOST_ASSERT( middle != end );

			//build 2 trees
			t.reset( new typename tree::element_type( ) );
			t->_min = *( boost::begin( sorted_points ) + begin );
			t->_max = *( boost::begin( sorted_points ) + end - 1 );
			operator()( //left
				t->_left,
				sorted_points,
				begin,
				middle );
			operator()( //right
				t->_right,
				sorted_points,
				middle,
				end );

			using namespace boost::lambda;
			using boost::lambda::_1;
			apply_to_sub_tree< d, traits, k >()(
				t,
				boost::lambda::bind< void >(
					create_tree< d, traits, k - 1 >(),
					boost::lambda::_1,
					boost::make_iterator_range(
						boost::begin( sorted_points ) + begin,
						boost::begin( sorted_points ) + end ) ) );
		}
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct create_tree< d, traits, 0 >
{

	template<
		typename Dummy,
		typename PointsRange
	>
	void operator()(
		Dummy & t,
		const PointsRange & points ) const
	{
	}
};

/** Checks whether points are in the correct child of a node. */
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k,
	unsigned i
>
struct bounds_checker
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef std::pair< point, weight > heaviest_point;

	point _min;
	point _max;
	bounds_checker( point min, point max ) : _min( min ) , _max( max ) { }
	void operator()( typename node< d, traits, i >::ptr & t ) const
	{
		if( t->_point )
		{
			typename traits::point_less_k cmp( k - 1 );
			BOOST_ASSERT( ! cmp( t->_point, _min ) );
			BOOST_ASSERT( ! cmp( _max, t->_point ) );
		}
		else
		{
			using namespace boost::lambda;
			using boost::lambda::_1;
			apply_to_sub_tree< d, traits, i >()(
				t,
				boost::lambda::bind< void >(
					bounds_checker< d, traits, k, i - 1 >( _min, _max ),
					boost::lambda::_1 ) );
		}
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,				/**< The traits of the problem. */
	unsigned k
>
struct bounds_checker< d, traits, k, 1 >
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef std::pair< point, weight > heaviest_point;

	point _min;
	point _max;
	bounds_checker( point min, point max ) : _min( min ) , _max( max ) { }
	void operator()( typename node< d, traits, 1 >::ptr & t ) const
	{
		if( t->_point )
		{
			typename traits::point_less_k cmp( k - 1 );
			BOOST_ASSERT( ! cmp( t->_point, _min ) );
			BOOST_ASSERT( ! cmp( _max, t->_point ) );
		}
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
struct check_tree;

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned m
>
struct check_sub_tree
{
	void operator()( typename node< d, traits, m >::ptr & t ) const
	{
		//must have sub tree xor be leaf node
		BOOST_ASSERT( ( bool( t->_point ) ) ^ ( bool( t->_sub_tree ) ) );

		//max weight must be the same as sub-tree's
		BOOST_ASSERT( t->_max_weight == t->_sub_tree->_max_weight );

		//check sub-tree has exactly same points as children do
		{
			typedef std::set< typename traits::point > point_set;

			point_set sub_tree_points;
			visit_tree< d, traits, m - 1 >()( t->_sub_tree, make_get_points< d, traits, m - 1 >( std::inserter( sub_tree_points, sub_tree_points.begin() ) ) );

			point_set tree_points;
			visit_tree< d, traits, m >()( t, make_get_points< d, traits, m >( std::inserter( tree_points, tree_points.begin() ) ) );

			BOOST_ASSERT( sub_tree_points == tree_points );
		}

		//recurse down the sub tree
		using namespace boost::lambda;
		using boost::lambda::_1;
		apply_to_sub_tree< d, traits, m >()(
			t,
			boost::lambda::bind< void >(
				check_tree< d, traits, m - 1 >(),
				boost::lambda::_1 ) );
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct check_sub_tree< d, traits, 1 >
{
	void operator()( typename node< d, traits, 1 >::ptr & t ) const
	{
		//nothing to do
	}
};


/** Checks trees are well formed. */
template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits,			/**< The traits of the problem. */
	unsigned k
>
struct check_tree
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef std::pair< point, weight > heaviest_point;

	void operator()( typename node< d, traits, k >::ptr & t ) const
	{
		if( t )
		{
			//if has left child, must have right and vice versa
			BOOST_ASSERT( ( bool( t->_left ) ) ^ ( bool( ! t->_right ) ) );

			//must be internal node xor be leaf node
			BOOST_ASSERT( ( bool( t->_point ) ) ^ ( bool( t->_left ) ) );

			//the _max can't be less than the _min
			BOOST_ASSERT( ! typename traits::point_less_k( k - 1 )( t->_max, t->_min ) );

			//is it a leaf node?
			if( ! t->_left )
			{
				//yes - _min must be the same as _max
				BOOST_ASSERT( t->_min == t->_max );
			}
			else
			{
				//check the _min _max bounds for all children
				visit_tree< d, traits, k >()( t->_left, bounds_checker< d, traits, k, k >( t->_min, t->_max ) );
				visit_tree< d, traits, k >()( t->_right, bounds_checker< d, traits, k, k >( t->_min, t->_max ) );

				//the left _max must be less than the right _min in the k'th dimension
				BOOST_ASSERT( typename traits::point_less_k( k - 1 )( t->_left->_max, t->_right->_min ) );

				//the _max_weight must be >= than the _max_weight of both children
				BOOST_ASSERT( t->_max_weight >= t->_left->_max_weight );
				BOOST_ASSERT( t->_max_weight >= t->_right->_max_weight );

				//check sub_tree properties
				check_sub_tree< d, traits, k >()( t );

				//recurse down both children
				operator()( t->_left );
				operator()( t->_right );
			}
		}
	}
};

template<
	unsigned d,					/**< Dimension of the range tree. */
	typename traits				/**< The traits of the problem. */
>
struct RangeTree
{
	typedef typename traits::weight weight;
	typedef typename traits::point point;
	typedef typename traits::point_less point_less;
	typedef typename traits::point_less_k point_less_k;
	typedef typename traits::point_equal point_equal;
	typedef typename traits::point_dominate point_dominate;
	typedef std::pair< point, weight > heaviest_point;


};

} // namespace range_tree_ns

using range_tree_ns::RangeTree;




namespace max_chain_ns {


template<
	unsigned d,				/**< Dimensions of the problem. */
	typename BoxT,				/**< Box type. */
	typename BoxTraits,			/**< Methods/typedefs to access box properties. */
	unsigned k
>
struct build_tree
{
	typedef BoxT box;
	typedef BoxTraits box_traits;
	typedef typename box_traits::box_ptr box_ptr;
	typedef typename box_traits::weight weight;
	typedef typename box_traits::coord coord;
	typedef typename box_traits::point point;
	typedef typename box_traits::point_equal point_equal;
	typedef typename box_traits::point_less point_less;
	typedef typename box_traits::point_less_k point_less_k;
	typedef typename box_traits::point_dominate point_dominate;
	typedef typename box_traits::data data;

	typedef std::vector< point > point_vec;
	typedef std::vector< point_vec > point_array;
	typedef std::set< point, point_less > point_set_equality;		/**< Order/equality by object contents. */
	typedef std::set< point > point_set;							/**< Order/equality by object address. */

	typedef std::vector< unsigned > idx_vec;
	typedef std::vector< idx_vec > idx_array;

	typedef std::multimap< point, box_ptr, point_less > point_to_box_multimap;
	typedef std::map< point, box_ptr, point_less > point_to_box_map;
	typedef std::pair< typename point_to_box_multimap::iterator, typename point_to_box_multimap::iterator > point_to_box_map_range;
	typedef std::map< box_ptr, weight > box_to_weight_map;
	typedef std::map< box_ptr, box_ptr > box_to_box_map;

	typedef RangeTree< d - 1, box_traits > range_tree;
	typedef typename range_tree_ns::create_tree< d, box_traits, d - 1 >::tree tree;

	point_set_equality points;
	point_set_equality end_points;
	tree t;
	point_to_box_multimap start_points_to_boxes;
	point_to_box_multimap end_points_to_boxes;
	box_to_weight_map chain_weights;					/**< The weight of the maximal chain that ends in this box. */
	point_to_box_map end_point_to_maximal_box;			/**< Points to the box that ends the maximal chain to this point. */
	typename box_traits::box_ptr_list & maximal_chain;

	build_tree( typename box_traits::box_ptr_list & mc ) : maximal_chain( mc ) { }

	template<
		typename BoxRange
	>
	weight operator()( const BoxRange & boxes )
	{

		/**
		Generate a set containing all the points in the problem.
		Some point objects will represent the same point. We only want to deal with one point
		object per actual point. We generate a map to deal with this.
		Also we need a map from points to boxes which start with that point and another map
		to boxes which end with that point.
		*/
		BOOST_FOREACH( const box & b, boxes )
		{
			typedef std::pair< typename point_set_equality::iterator, bool > insert_result;
			box_ptr box_p = boost::addressof( b );

			//start
			{
				point p = box_traits::get_start( b );
				points.insert( p );
				start_points_to_boxes.insert( typename point_to_box_multimap::value_type( p, box_p ) );
			}

			//end
			{
				point p = box_traits::get_end( b );
				points.insert( p );
				end_points.insert( p );
				end_points_to_boxes.insert( typename point_to_box_multimap::value_type( p, box_p ) );
			}
		}

		range_tree_ns::template create_tree< d, box_traits, d - 1 >()( t, end_points );
		//range_tree_ns::visit_tree< d, box_traits, d - 1 >()( t, range_tree_ns::print_point< d, box_traits, d - 1 >() );
		range_tree_ns::template visit_tree< d, box_traits, d - 1 >()(
			t,
			range_tree_ns::template check_tree< d, box_traits, d - 1 >() );

		range_tree_ns::template insert_point< d, box_traits, d - 1 > point_inserter;

		//for each point (ordered in the d'th dimension)
		coord last_dth_coord = std::numeric_limits< coord >::min();
		BOOST_FOREACH( point p, points )
		{
			//range_tree_ns::print_point< d >().print( p );

			//check points actually are ordered in d'th dimension...
			if( last_dth_coord > ( *p )[ d - 1 ] )
			{
				throw std::logic_error( "box_traits::point_less must order points based on d'th dimension" );
			}
			last_dth_coord = ( *p )[ d - 1 ];

			//get the heaviest point dominated by this point
			typename range_tree::heaviest_point max_weight_dominated_by_p = get_heaviest_point( p );

			//is this point an end-point?
			point_to_box_map_range end_boxes = end_points_to_boxes.equal_range( p );
			if( end_boxes.first != end_boxes.second )
			{
				//our best box so far for this end point is the box that ends at the previous
				//heaviest point
				box_ptr prev_box = get_maximal_box_for( max_weight_dominated_by_p.first );
				point_inserter( t, p, max_weight_dominated_by_p.second );
				if( 0 != prev_box )
				{
					end_point_to_maximal_box[ p ] = prev_box;
				}

				//for all boxes, b, that end at this point, p, see if we can find a better chain...
				BOOST_FOREACH(
					typename point_to_box_multimap::value_type v,
					boost::make_iterator_range(
						end_boxes.first,
						end_boxes.second ) )
				{
					box_ptr b = v.second;

					//get the weight of the maximal chain that ends in this box
					typename box_to_weight_map::const_iterator w_b = chain_weights.find( b );
					BOOST_ASSERT( chain_weights.end() != w_b ); //we must have filled this in previously
					const weight box_maximal_weight = w_b->second;

					//is the heaviest point lighter than the box and any other boxes we have inserted here?
					if( max_weight_dominated_by_p.second < box_maximal_weight )
					{
						//yes so update this point as having the weight of this box's maximal chain
						point_inserter( t, p, box_maximal_weight );

						//update our record of the heaviest box for this point
						max_weight_dominated_by_p.second = box_maximal_weight;
						max_weight_dominated_by_p.first = p;

						//update our map that says that the best box that ends here is this box
						end_point_to_maximal_box[ p ] = b;

#ifdef _DEBUG
						BOOST_ASSERT( ( typename range_tree_ns::max_weight< d, box_traits, d - 1 >()( t, p ).second == box_maximal_weight ) );
						//check the weight of the maximal chain to the point agrees with the weight in the tree.
						check_maximal_chain_weight( max_weight_dominated_by_p.second, max_weight_dominated_by_p.first );
#endif //_DEBUG
					}
				}
			}

#ifdef _DEBUG
			//check the weight of the maximal chain to the point agrees with the weight in the tree.
			check_maximal_chain_weight( max_weight_dominated_by_p.second, max_weight_dominated_by_p.first );
#endif //_DEBUG

			//for all boxes, b, that start at this point, p.
			point_to_box_map_range start_boxes = start_points_to_boxes.equal_range( p );
			BOOST_FOREACH(
				typename point_to_box_multimap::value_type v,
				boost::make_iterator_range(
					start_boxes.first,
					start_boxes.second ) )
			{
				box_ptr b = v.second;

				//the weight of the maximal chain ending in this box is the weight of the box
				//plus the weight of the heaviest point before it..
				const weight chain_weight = max_weight_dominated_by_p.second + b->_weight;

				//store the maximal chain weight
				chain_weights.insert( typename box_to_weight_map::value_type( b, chain_weight ) );
			}

#ifdef _DEBUG
			//check the weight of the maximal chain to the point agrees with the weight in the tree.
			check_maximal_chain_weight( max_weight_dominated_by_p.second, max_weight_dominated_by_p.first );
#endif //_DEBUG

		}

		//find the heaviest point in the range tree
		const typename range_tree::heaviest_point overall_max = get_heaviest_point( box_traits::always_dominates_point() );
		//std::cout << "Max weight = " << overall_max.second << "\n";

		//did we find a point?
		if( ! point_equal()( box_traits::always_dominated_point(), overall_max.first ) )
		{
			//find the box that ends at this point
			box_ptr maximal_box = get_maximal_box_for( overall_max.first );

			//get the maximal chain to this point
			get_maximal_chain_to(
				maximal_box,
				std::front_inserter( maximal_chain ) );
		}

		return overall_max.second;
	}

	//check the weight we have for the maximal chain actually is correct...
	void
	check_maximal_chain_weight( weight maximal_chain_weight, point p )
	{
		//this is recursive - so if already checking just ignore
		static bool already_checking = false;
		if( already_checking )
		{
			return;
		}
		already_checking = true;

		if( 0.0 != maximal_chain_weight )
		{
			//did we have a point?
			if( point_equal()( box_traits::always_dominated_point(), p ) )
			{
				//no
				BOOST_ASSERT( weight( 0.0 ) == maximal_chain_weight );
			}
			else
			{
				//we found a point
				typename box_traits::box_ptr_list max_chain_to_point;
				typename point_to_box_map::const_iterator b_it = end_point_to_maximal_box.find( p );
				if( end_point_to_maximal_box.end() == b_it )
				{
					BOOST_ASSERT( end_point_to_maximal_box.end() != b_it );
				}
				get_maximal_chain_to(
					b_it->second,
					std::front_inserter( max_chain_to_point ) );
				weight w = weight( 0.0 );
				BOOST_FOREACH( box_ptr b, max_chain_to_point )
				{
					w += b->_weight;
				}
				const bool is_close =
					boost::test_tools::check_is_close(
						w,
						maximal_chain_weight,
						BIO_FPC_NS::percent_tolerance( 0.01 )
					);
				BOOST_ASSERT( is_close );
			}
		}
		already_checking = false;
	}

	/** Gets the maximal chain to the given box. */
	template<
		typename FrontInsertIt
	>
	void
	get_maximal_chain_to(
		box_ptr b,
		FrontInsertIt output_it )
	{
		//use the maximal chain previous box map to reverse down the maximal chain
		//keep going until we cannot find any more previous boxes
		while( 0 != b )
		{
			//output this box
			*output_it = b;
			++output_it;

			//go to the previous
			b = get_previous_box_for( box_traits::get_start( *b ) );
		}
	}

	/** Gets the heaviest box that is dominated by p. */
	box_ptr
	get_previous_box_for(
		point p )
	{
		return get_maximal_box_for( get_heaviest_point( p ).first );
	}

	/** Get the heaviest point that is dominated by this point. */
	typename range_tree::heaviest_point
	get_heaviest_point( point p )
	{
		//get the heaviest point dominated by this point
		const typename range_tree::heaviest_point heavy_point = range_tree_ns::template max_weight< d, box_traits, d - 1 >()( t, p );

#ifdef _DEBUG
		//check the weight of the maximal chain from the point agrees with the weight in the tree.
		check_maximal_chain_weight( heavy_point.second, heavy_point.first );
#endif

		return heavy_point;
	}

	/** Gets the box that ends at p that is in the maximal chain to p. I.e. the best box. */
	box_ptr
	get_maximal_box_for(
		point p )
	{
		typename point_to_box_map::const_iterator b_it = end_point_to_maximal_box.find( p );
		return
			b_it != end_point_to_maximal_box.end()
				? b_it->second
				: 0;
	}

};

template<
	unsigned d,				/**< Dimensions of the problem. */
	typename BoxT,				/**< Box type. */
	typename BoxTraits			/**< Methods/typedefs to access box properties. */
>
struct build_tree< d, BoxT, BoxTraits, 1 >
{
	typedef typename BoxTraits::weight weight;

	typename BoxTraits::box_ptr_list & maximal_chain;

	build_tree( typename BoxTraits::box_ptr_list & mc ) : maximal_chain( mc ) { }

	template<
		typename BoxRange
	>
	weight operator()( const BoxRange & boxes )
	{
		throw std::logic_error( "MaxChain algorithm not implemented for single sequences" );

		return weight( 0.0 );
	}
};

template<
	unsigned d,				/**< Dimensions of the problem. */
	typename BoxT,				/**< Box type. */
	typename BoxTraits			/**< Methods/typedefs to access box properties. */
>
struct build_tree< d, BoxT, BoxTraits, 0 >
{
	typedef typename BoxTraits::weight weight;

	typename BoxTraits::box_ptr_list & maximal_chain;

	build_tree( typename BoxTraits::box_ptr_list & mc ) : maximal_chain( mc ) { }

	template<
		typename BoxRange
	>
	weight operator()( const BoxRange & boxes )
	{
		return weight( 0.0 );
	}
};


template<
	unsigned d,				/**< Dimensions of the problem. */
	typename BoxT,				/**< Box type. */
	typename BoxTraits			/**< Methods/typedefs to access box properties. */
>
struct MaxChain
{
	typedef BoxT box;
	typedef BoxTraits box_traits;
	typedef typename box_traits::box_ptr box_ptr;
	typedef typename box_traits::weight weight;
	typedef typename box_traits::coord coord;
	typedef typename box_traits::point point;
	typedef typename box_traits::point_equal point_equal;
	typedef typename box_traits::point_less point_less;
	typedef typename box_traits::point_less_k point_less_k;
	typedef typename box_traits::point_dominate point_dominate;
	typedef typename box_traits::data data;

	typedef std::vector< point > point_vec;
	typedef std::vector< point_vec > point_array;
	typedef std::set< point, point_less > point_set_equality;		/**< Order/equality by object contents. */
	typedef std::set< point > point_set;							/**< Order/equality by object address. */

	typedef std::vector< unsigned > idx_vec;
	typedef std::vector< idx_vec > idx_array;

	typedef std::multimap< point, box_ptr, point_less > point_to_box_multimap;
	typedef std::map< point, box_ptr, point_less > point_to_box_map;
	typedef std::pair< typename point_to_box_multimap::iterator, typename point_to_box_multimap::iterator > point_to_box_map_range;
	typedef std::map< box_ptr, weight > box_to_weight_map;
	typedef std::map< box_ptr, box_ptr > box_to_box_map;

	typedef RangeTree< d - 1, box_traits > range_tree;

	typename box_traits::box_ptr_list maximal_chain;
	weight maximal_weight;


	template<
		typename BoxRange
	>
	MaxChain( const BoxRange & boxes )
	: maximal_weight( 0.0 )
	{
		if( box_traits::get_dimensions() != d )
		{
			throw std::invalid_argument( "Wrong number of dimensions" );
		}

		build_tree< d, box, box_traits, d > tree_builder( maximal_chain );
		maximal_weight = tree_builder( boxes );
	}

	template< unsigned k >
	struct calculate_point_indices
	{
		template<
			typename PointsArray,
			typename IndexArray
		>
		void operator()(
			const PointsArray & sorted_points,
			IndexArray & point_indices ) const
		{
			point_indices.resize( k - 1 );
			for( unsigned i = 1; k > i; ++i )
			{
				BOOST_FOREACH( point p, sorted_points[ i ] )
				{
					typedef typename point_vec::const_iterator it;

					//find the point in the dimension before
					point_vec & prev_dim( sorted_points[ i - 1 ] );
					it lower =
						std::lower_bound(
							prev_dim.begin(),
							prev_dim.end(),
							p,
							point_less_k( i - 1 ) );

					it upper =
						std::upper_bound(
							prev_dim.begin(),
							prev_dim.end(),
							p,
							point_less_k( i - 1 ) );

					it p_it = std::find_if( lower, upper, boost::bind< bool >( point_equal(), p, _1 ) );

					const unsigned idx = p_it - prev_dim.begin();
					point_indices[ i - 1 ].push_back( idx );
				}
			}
		}
	};

};


/**
Generic box type for use with max chain algorithm.
*/
template<
	unsigned d,
	typename Data,
	typename Coord,
	typename Weight
>
struct MaxChainBox
{
	struct traits
		: point_traits< Coord, d >
	{
		typedef typename point_traits< Coord, d >::point point;
		typedef typename point_traits< Coord, d >::coord_array coord_array;
		typedef typename point_traits< Coord, d >::coord coord;

		typedef Weight weight;
		typedef Data data;
		typedef MaxChainBox< d, data, coord, weight > box;
		typedef const box * box_ptr;
		typedef std::vector< box > box_vec;					/**< Vector of boxes. */
		typedef std::list< box_ptr > box_ptr_list;			/**< List of box pointers. */

		//typedef MaxChainBox< d, coord, data > box;
		//typedef typename traits::box_ptr box_ptr;
		//typedef typename point_traits< Coord, d >::coord_array coord_array;
		//typedef typename point_traits< Coord, d >::point point;

		static unsigned get_dimensions() { return d; }

		static box make_box( data _d, weight w, coord_array s, coord_array e ) { return box( _d, w, s, e ); }
		static weight get_weight( const box & b ) { return b._weight; }
		static data get_data( const box & b ) { return b._data; }
		static point get_start( const box & b ) { return boost::addressof( b._start ); }
		static point get_end( const box & b ) { return boost::addressof( b._end ); }
	};


	MaxChainBox(
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





template<
	unsigned d,
	typename value_traits
>
struct max_chain_algorithm
{
	typedef MaxChainBox<
		d,
		typename value_traits::data,
		typename value_traits::coord,
		typename value_traits::weight
	> box;


	typedef MaxChain<
		d,
		box,
		typename box::traits
	> algorithm;



	template<
		typename ValueSequenceRange,
		typename OutputIt
	>
	bool
	operator()(
		const ValueSequenceRange & hits,
		OutputIt output_it,
		unsigned box_limit = 0 ) const
	{
	    //
	    // Calculate the boxes from the values
	    //
	    typedef typename box::traits box_traits;
        typedef typename box_traits::box_vec box_vec;
        box_vec boxes;
        auto values = box_generator< value_traits >::template boxes_from_sequences< box_traits >(
            hits,
            std::back_inserter( boxes ),
            box_limit
        );
		log_stream() << "Calculating max chain: # boxes = " << boxes.size() << std::endl;

		bool ran_algorithm = false;

		//do we have a limit or do we have fewer boxes than it?
		if( ! box_limit || boxes.size() < box_limit )
		{
			algorithm alg( boxes );

			BOOST_FOREACH( typename algorithm::box_ptr b, alg.maximal_chain )
			{
				*output_it = b;
				++output_it;
			}
			ran_algorithm = true;
		}

		return ran_algorithm;
	}
};



struct print_box
{
	std::ostream & _os;
	print_box( std::ostream & os ) : _os( os ) { }
	template< typename BoxPtr >
	void operator()( BoxPtr b )
	{
		std::copy( b->_start.begin(), b->_start.end(), std::ostream_iterator< int >( _os, "," ) );
		_os << "; ";
		std::copy( b->_end.begin(), b->_end.end(), std::ostream_iterator< int >( _os, "," ) );
		_os << "; " << b->_weight << "; " << *( b->_data ) << std::endl;
	}
};



/**
A macro for each branch of the switch on the problem dimension.
*/
#define BOOST_PP_LOCAL_MACRO(n) \
	case n: \
		return \
			max_chain_algorithm< n, ValueTraits >()(  \
				hits,  \
				output_it, \
				num_boxes_limit );

#define BOOST_PP_LOCAL_LIMITS (2, BIO_MAX_CHAIN_MAX_SEQUENCES)


/**
Calculates the maximal chain that is a subsequence of each sequence of hits in the HitSequenceRange.
*/
template<
	typename ValueTraits,
	typename ValueSequenceRange,
	typename OutputIt
>
bool
max_chain(
	const ValueSequenceRange & hits,
	OutputIt output_it,
	unsigned num_boxes_limit = 0 )
{
	using namespace boost::lambda;
	using boost::lambda::_1;

	switch( boost::size( hits ) )
	{
	case 0:
		return false; //nothing to do if we don't have any sequences

	case 1: //cannot handle single sequences yet
		throw
			std::logic_error(
				BIO_MAKE_STRING(
					"Max chain algorithm not implemented for " << boost::size( hits ) << " sequences. Limit is " << BIO_MAX_CHAIN_MAX_SEQUENCES ) );

#include BOOST_PP_LOCAL_ITERATE() //expands to fill in cases for range defined above

	default: //too many sequences to handle
		throw
			std::logic_error(
				BIO_MAKE_STRING(
					"Max chain algorithm not implemented for " << boost::size( hits ) << " sequences. Limit is " << BIO_MAX_CHAIN_MAX_SEQUENCES ) );
	}
}

#undef BIO_max_chain_case_instance




} // namespace max_chain_ns

using max_chain_ns::max_chain;
using max_chain_ns::MaxChain;


BIO_NS_END

#endif //BIO_MAX_CHAIN_H_

