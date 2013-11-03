/**
@file

Copyright John Reid 2013
*/

#ifndef BIO_MAX_CHAIN_BOXES_H_
#define BIO_MAX_CHAIN_BOXES_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER


#include <bio/defs.h>

#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/foreach.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm_ext/is_sorted.hpp>
#include <boost/range/adaptor/indirected.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>


#ifndef BIO_MAX_CHAIN_MAX_SEQUENCES
# define BIO_MAX_CHAIN_MAX_SEQUENCES 5
#endif //BIO_MAX_CHAIN_MAX_SEQUENCES


BIO_NS_START



/**
 * Iterates across multiple dimensions. Used to generate combinations of values to make up boxes.
 */
struct multi_iterator
{
    typedef std::vector< unsigned > index;

    /** Construct from range of ranges. (Need to be model of ForwardRange). */
    template<
        typename RangeOfRanges
    >
    multi_iterator( const RangeOfRanges & range_of_ranges )
    : _ind( boost::size( range_of_ranges ), 0 )
    , _begin( boost::size( range_of_ranges ), 0 )
    {
        BOOST_FOREACH( const typename boost::range_value< RangeOfRanges >::type & range, range_of_ranges )
        {
            _end.push_back( boost::size( range ) );
        }
    }

    /** Construct from a beginning and an end index. */
    multi_iterator( const index & begin, const index & end )
        : _ind( begin )
        , _begin( begin )
        , _end( end )
    {
        if( begin.size() != end.size() )
        {
            throw std::invalid_argument( "begin and end sizes are different" );
        }
    }
    index _ind;
    index _begin;
    index _end;

    unsigned get_dimensions() const
    {
        return _ind.size();
    }

    bool at_end() const
    {
        for( unsigned d = 0; get_dimensions() != d; ++d )
        {
            if( _ind[ d ] == _end[ d ] )
            {
                return true;
            }
        }
        return false;
    }

    void next()
    {
        for( unsigned d = 0; get_dimensions() != d; ++d )
        {
            _ind[ d ]++;
            if( _end[ d ] == _ind[ d ] )
            {
                _ind[ d ] = _begin[ d ];
            }
            else
            {
                return; //we have incremented the index
            }
        }
        //set ind to be the end
        _ind = _end;
    }
};



template<
    typename Coord,        /**< The coordinate type of the points. */
    unsigned d              /**< The dimensions of the points. */
>
struct point_traits
{
    typedef Coord coord;
    typedef boost::array< coord, d > coord_array;
    typedef const coord_array * point;

    struct point_equal
    {
        bool operator()( point lhs, point rhs ) const
        {
            for( unsigned k = 0; d != k; ++k )
            {
                if( ( *lhs )[ k ] != ( *rhs )[ k ] )
                {
                    return false;
                }
            }
            return true;
        }
    };

    /** Compare 2 points one dimension at a time. Equal points dominate themselves. */
    struct point_dominate
    {
        bool operator()( point dominatrix, point slave ) const
        {
            for( unsigned k = 0; d != k; ++k )
            {
                if( ( *dominatrix )[ k ] < ( *slave )[ k ] )
                {
                    return false;
                }
            }
            return true;
        }
    };

    /** Compare 2 points one dimension at a time. Use d'th dimension as first to check. */
    struct point_less
    {
        bool operator()( point lhs, point rhs ) const
        {
            for( unsigned k = d; 0 != k; --k )
            {
                if( ( *lhs )[ k - 1 ] == ( *rhs )[ k - 1 ] )
                {
                    continue;
                }
                else
                {
                    return ( *lhs )[ k - 1 ] < ( *rhs )[ k - 1 ];
                }
            }
            return false; //they are equal
        }
    };

    /** Compares 2 points just in dimension k, then k+1%d, then k+2%d, etc... */
    struct point_less_k
    {
        unsigned _k;
        point_less_k( unsigned k ) : _k( k ) { }
        bool operator()( point lhs, point rhs ) const
        {
            for( unsigned i = 0; d != i; ++i )
            {
                const unsigned j = ( i + _k ) % d;
                if( ( *lhs )[ j ] == ( *rhs )[ j ] )
                {
                    continue;
                }
                return ( *lhs )[ j ] < ( *rhs )[ j ];
            }
            return false; //are equal
        }
    };

    /**< A point that is dominated by all others. */
    static point always_dominated_point()
    {
        static boost::shared_ptr< coord_array > a;
        if( ! a )
        {
            a.reset( new coord_array );
            std::fill( a->begin(), a->end(), std::numeric_limits< coord >::min() );
        }
        return a.get();
    }

    /**< A point that dominates all others. */
    static point always_dominates_point()
    {
        static boost::shared_ptr< coord_array > a;
        if( ! a )
        {
            a.reset( new coord_array );
            std::fill( a->begin(), a->end(), std::numeric_limits< coord >::max() );
        }
        return a.get();
    }
};



/**
 * This code has been re-factored into the box_generator class below...
 */
//template<
//    typename ValueTraits,               /**< Traits to access members of values. */
//    typename BoxTraits,                 /**< Our output type traits. */
//    typename SequencesRange,           /**< A range of ranges of values. */
//    typename OutIt                      /**< Output iterator to write boxes to. */
//>
//void get_boxes_for(
//    const SequencesRange & sequences,
//    OutIt output_it )
//{
//    using namespace boost;
//
//    typedef ValueTraits value_traits;
//    typedef typename range_value< SequencesRange >::type sequence;
//    typedef typename range_value< sequence >::type value;
//    typedef typename value_traits::character character;
//    typedef std::set< character > char_set;
//
//    typedef OutIt output_iterator;
//
//    typedef BoxTraits box_traits;
//    typedef typename box_traits::weight weight;
//    typedef typename box_traits::data data;
//    typedef typename box_traits::coord_array point;
//
//    const unsigned dimensions = size( sequences );
//    if( point::static_size != dimensions )
//    {
//        throw std::invalid_argument(
//            BIO_MAKE_STRING( "Wrong # of input sequences: " << unsigned( point::static_size ) << " != " << dimensions ) );
//    }
//
//    if( empty( sequences ) )
//    {
//        return;
//    }
//
//    //for each value in first sequence.
//    char_set chars;
//    BOOST_FOREACH( const value & v, sequences[ 0 ] )
//    {
//        const character c = ValueTraits::get_char( v );
//
//        //have we looked at this char already
//        if( chars.end() != chars.find( c ) )
//        {
//            //yes - so ignore and continue...
//            continue;
//        } else {
//            chars.insert( c );
//        }
//
//        typedef std::vector< const value * > value_ptr_vec;
//        typedef std::vector< value_ptr_vec > value_ptr_array;
//        value_ptr_array values;
//        BOOST_FOREACH( const sequence & s, sequences ) //for each input sequence
//        {
//            values.push_back( value_ptr_vec() );
//            BOOST_FOREACH( const value & v, s ) //for each value in input sequence
//            {
//                if( ValueTraits::get_char( v ) == c )
//                {
//                    values.rbegin()->push_back( boost::addressof( v ) );
//                }
//            }
//        }
//
//        //now we have all the values of the given character, we must iterate across
//        //all possible combinations
//        for( multi_iterator it( values ); ! it.at_end(); it.next() )
//        {
//            //calculate the box
//            point start;
//            point end;
//            weight w = weight( 0.0 );
//
//            //for each dimension
//            for( unsigned d = 0; dimensions != d; ++d )
//            {
//                const value * v = values[ d ][ it._ind[ d ] ];
//                start[ d ] = value_traits::get_start( *v );
//                end[ d ] = value_traits::get_end( *v );
//                w += value_traits::get_weight( *v );
//            }
//            w /= dimensions;
//
//            *output_it =
//                box_traits::make_box(
//                    value_traits::get_data( *( values[ 0 ][ it._ind[ 0 ] ] ) ), //get the data from the first value...
//                    w,
//                    start,
//                    end );
//            ++output_it;
//        }
//    }
//}


/**
 * Calculate the product of a range.
 */
template< typename Range >
typename boost::range_value< Range >::type
product( const Range & r ) {
    typename boost::range_value< Range >::type result = 1;
    BOOST_FOREACH( typename boost::range_value< Range >::type x, r ) {
        result *= x;
    }
    return result;
}


template< typename ValueTraits >
struct box_generator {
    typedef typename ValueTraits::hit value;
    typedef typename ValueTraits::weight weight;
    typedef typename ValueTraits::character character;
    typedef std::vector< value > value_vec;
    typedef std::vector< value_vec > value_vec_vec;
    typedef boost::shared_ptr< value_vec_vec > value_vec_vec_ptr;
    typedef std::map< character, value_vec_vec_ptr > map_by_char;
    typedef typename map_by_char::iterator map_by_char_iterator;
    typedef boost::iterator_range< typename boost::range_iterator< const value_vec >::type > thresholded_values;
//    typedef boost::transformed_range< thresholded_values > thresholded_value_vecs;




    /// Compares the weight of 2 values
    struct weight_less_than {
        typedef bool result_type;
        bool
        operator()( const value & lhs, const value & rhs ) const {
            return ValueTraits::get_weight( lhs ) < ValueTraits::get_weight( rhs );
        }
        bool
        operator()( weight w, const value & rhs ) const {
            return w < ValueTraits::get_weight( rhs );
        }
    };


    /// Is the weight of a value higher than a threshold?
    struct threshold_weight_predicate {
        typedef bool result_type;
        bool
        operator()( weight threshold, const value & v ) const {
            return ValueTraits::get_weight( v ) > threshold;
        }
    };


    /**
     * Return a range that contains only those values above the threshold.
     */
    static
    thresholded_values
    threshold_values( const value_vec & values, weight threshold ) {
        return boost::make_iterator_range(
            boost::upper_bound( values, threshold, threshold_weight_predicate() ),
            values.end()
        );
    }


//    /**
//     * Return a range that contains only those values above the threshold.
//     */
//    static
//    thresholded_value_vecs
//    threshold_value_vecs( const value_vec_vec & values, weight threshold ) {
//        using boost::lambda::bind;
//        using boost::lambda::_1;
//        using boost::adaptors::transformed;
//        return values | transformed( bind( &threshold_values, _1, threshold ) );
//    }


    /**
     * Calculate the number of boxes the values will generate at the given threshold.
     * This is the product of the number of values above the threshold in each value vector.
     */
    static
    unsigned
    calculate_num_boxes( const value_vec_vec & values, weight threshold ) {
        using boost::adaptors::transformed;
        using boost::lambda::_1;
        return product(
            values
            | transformed( bind( &threshold_values, _1, threshold ) )
            | transformed( boost::size< thresholded_values >( _1 ) )
        );
    }


    /**
     * Take the values and create boxes from them.
     * The boxes contain pointers to the values returned in the map_by_char
     * container so this container's lifetime should be maintained for as long
     * as the boxes lifetimes.
     */
    template<
        typename BoxTraits,       /**< Our output type traits. */
        typename SequencesRange, /**< A range of ranges of values. */
        typename OutIt            /**< Output iterator to write boxes to. */
    >
    static
    boost::shared_ptr< map_by_char >
    boxes_from_sequences(
        const SequencesRange & sequences,
        OutIt output_it,
        unsigned max_num_boxes = 0
    ) {
        boost::shared_ptr< map_by_char > values_by_char( new map_by_char );
        box_generator< ValueTraits >::template organise_values_by_character( sequences, *values_by_char );
        if( max_num_boxes && ! values_by_char->empty() ) { // if we have a limit and some values
            num_boxes_bounder( *values_by_char, max_num_boxes )();
        }
        box_generator< ValueTraits >::template get_boxes_for< BoxTraits >( *values_by_char, output_it );
        return values_by_char;
    }


    /**
     * Take the values and separate them into a data structure where they are organised by character.
     */
    template<
        typename SequencesRange             /**< A range of ranges of values. */
    >
    static
    void
    organise_values_by_character(
        const SequencesRange & sequences,
        map_by_char & result
    ) {
        using namespace boost;

        //
        // Examine value type of sequences
        //
        typedef typename range_value< SequencesRange >::type sequence;
        typedef typename range_value< sequence >::type seq_value;
        typedef is_same< value, seq_value > are_values_same_type;
        BOOST_STATIC_ASSERT_MSG( are_values_same_type::value, "Value types should be equal." );

        //
        // Empty the result collection
        //
        result.clear();

        //
        // for each sequence
        //
        unsigned s = 0;
        for(
            typename range_iterator< const SequencesRange >::type seq_it = boost::begin( sequences );
            boost::end( sequences ) != seq_it;
            ++seq_it, ++s
        ) {

            //
            // for each value in the sequence
            //
            BOOST_FOREACH( const value & v, *seq_it ) {

                //
                // which character is the value for?
                //
                const character c = ValueTraits::get_char( v );

                //
                // Check if this character is in the results map
                //
                map_by_char_iterator char_it = result.find( c );
                if( char_it == result.end() ) {
                    //
                    // It isn't. If this is the first sequence we should
                    // initialise the data structures. If not then we
                    // will ignore it as we are only interested if there is an entry
                    // for the character already. I.e. it was in all
                    // the previous sequences.
                    //
                    if( 0 == s ) {
                        //
                        // This is the first sequence, initialise
                        //
                        char_it = result.insert(
                            typename map_by_char::value_type(
                                c,
                                value_vec_vec_ptr( new value_vec_vec )
                            )
                        ).first;
                        char_it->second->resize( size( sequences ) );
                    } else {
                        //
                        // This is not the first sequence, and there is no entry for the
                        // character in the results so we are not interested in it
                        //
                        continue;
                    }
                }

                //
                // add the value to the character's sequence
                //
                ( *char_it->second )[ s ].push_back( v );
            }

            //
            // We check which characters actually had values in this sequence as
            // we are not interested in any characters that did not have values.
            // We can remove them
            // from the result. Also we want to sort the vectors of values
            // to speed processing later.
            //
            std::vector< map_by_char_iterator > to_remove;
            for(
                map_by_char_iterator char_it = result.begin();
                result.end() != char_it;
                ++char_it
            ) {
                if( empty( ( *char_it->second )[ s ] ) ) {
                    to_remove.push_back( char_it ); // there were no values for this character so remove it from results
                } else {
                    sort( ( *char_it->second )[ s ], weight_less_than() ); // sort by weight
                }
            }
            BOOST_FOREACH( map_by_char_iterator char_it, to_remove ) {
                result.erase( char_it );
            }
        }
    }


    /**
     * Bounds the number of boxes by finding the correct threshold.
     */
    struct num_boxes_bounder {
        typedef typename value_vec::const_iterator iterator; ///< An iterator into the most nested collection of values, i.e. those values for one character in one sequence.
        struct sequence {
            value_vec * values;
            iterator lower, upper;
            sequence( value_vec & values )
            : values( &values )
            , lower( values.begin() )
            , upper( values.end() )
            {
                BOOST_ASSERT( boost::is_sorted( values, weight_less_than() ) );
            }
            auto make_range() const -> decltype( boost::make_iterator_range( lower, upper ) ) {
                return boost::make_iterator_range( lower, upper );
            }
            unsigned upper_count() const { return values->end() - upper; }
            unsigned lower_count() const { return values->end() - lower; }
        };
        typedef std::vector< sequence > sequence_vec; ///< A container of sequences.
        typedef std::map< character, sequence_vec > seqs_by_char; ///< A map from characters into containers of sequences.
        typedef std::vector< iterator > iterator_vec; ///< A pair of iterators representing an upper and lower bound.
        typedef std::map< character, iterator_vec > iterators_by_char; ///< A map from characters into containers of iterators.
        typedef boost::tuple< unsigned, character, unsigned, unsigned > char_count; ///< Counts for a character. tuple := (# boxes allowed, char, # boxes at upper, # boxes at lower)
        typedef std::vector< char_count > char_count_vec;
        typedef typename char_count_vec::iterator char_count_it;

        struct can_remove {
            typedef unsigned result_type;
            unsigned
            operator()( const char_count & count ) const {
                return count.template get< 0 >() - count.template get< 2 >();
            }
        };

        struct counts_less_than {
            typedef bool result_type;
            bool
            operator()( const char_count & lhs, const char_count & rhs ) const {
                return can_remove()( lhs ) < can_remove()( rhs );
            }
        };

        template< typename T >
        struct get_size {
            typedef unsigned result_type;
            unsigned
            operator()( const T & x ) const {
                //return x.size();
                return boost::size( x );
            }
        };

        struct get_count {
            typedef unsigned result_type;
            unsigned &
            operator()( char_count & x ) const {
                return x.template get< 0 >();
            }
        };

        seqs_by_char seqs; ///< The sequences of original values indexed by character.
        unsigned max_num_boxes; ///< Limit on the number of boxes.
        weight upper_weight; ///< Current upper weight. Using this weight is guaranteed to create the right number or fewer boxes.
        weight lower_weight; ///< Current lower weight. Using this weight is guaranteed to create too many boxes.
        weight max_weight; ///< Maximum weight in original values.
        weight min_weight; ///< Minimum weight in original values.


        /**
         * Check (only when BOOST_ASSERT is on) that some invariants that should hold actually do.
         */
        void
        check_invariants() {
            BOOST_ASSERT( calculate_num_boxes( false ) >=  max_num_boxes );
            BOOST_ASSERT( calculate_num_boxes( true  ) <  max_num_boxes );
        }


        /**
         * Reduce the number of values such that we are under a given upper limit of boxes as the
         * maximal chain algorithm is slow (quadratic?) in the number of boxes. The input values should
         * be sorted by weight as organise_values_by_character() does.
         */
        void
        operator()() {
            //
            // check there is something to restrict
            //
            if( calculate_num_boxes( false ) > max_num_boxes ) {
                //
                // If all the values have the same weight, there is nothing to do.
                // I.e. we are happy with the upper and lower weights set in the
                // constructor.
                //
                if( max_weight != min_weight ) {
                    restrict_bounds();
                }
                check_invariants();
                reduce_values();
            }
        }


        num_boxes_bounder( map_by_char & values_by_char, unsigned max_num_boxes )
        : max_num_boxes( max_num_boxes )
        , max_weight( -std::numeric_limits< weight >::max() ) // initialise to lowest  possible value
        , min_weight(  std::numeric_limits< weight >::max() ) // initialise to highest possible value
        {
            BOOST_ASSERT( ! values_by_char.empty() );

            //
            // Initialise our containers of pointers into the original values.
            // Also find the largest weight and the smallest.
            //
            BOOST_FOREACH( typename map_by_char::value_type values, values_by_char ) {
                BOOST_FOREACH( value_vec & vs, *values.second ) {
                    seqs[ values.first ].push_back( sequence ( vs ) );
                    if( ! vs.empty() ) {
                        const weight low  = ValueTraits::get_weight( vs.front() );
                        const weight high = ValueTraits::get_weight( vs.back()  );
                        BOOST_ASSERT( low <= high ); // weak test that values are sorted by weight
                        min_weight = std::min( min_weight, low  );
                        max_weight = std::max( max_weight, high );
                    }
                }
            }

            //
            // adjust lower weight so is just lower than the smallest weight
            //
            const weight weight_range = max_weight - min_weight;
            BOOST_ASSERT( weight_range >= 0. );
            weight adjustment = weight_range / 100.;
            if( 0. == adjustment ) { adjustment = weight_range; } // in case of rounding errors
            if( 0. == adjustment ) { adjustment = 1.e-3; } // in case of no range
            lower_weight = min_weight - adjustment;
            upper_weight = max_weight + adjustment;
        }


        /**
         * Calculate how many boxes with current bound (upper or lower).
         */
        unsigned
        calculate_num_boxes( bool use_upper ) {
            unsigned total_num_boxes = 0;
            BOOST_FOREACH( typename seqs_by_char::value_type seqs_for_char, seqs ) {
                total_num_boxes += calculate_num_boxes( seqs_for_char.second, use_upper );
            }
            return total_num_boxes;
        }


        /**
         * Calculate how many boxes with current bound (upper or lower).
         */
        unsigned
        calculate_num_boxes( const sequence_vec & seqs_for_char, bool use_upper ) {
            unsigned num_boxes_for_char = 1;
            BOOST_FOREACH( const sequence & seq, seqs_for_char ) {
                num_boxes_for_char *= seq.values->end() - ( use_upper ? seq.upper : seq.lower );
            }
            return num_boxes_for_char;
        }


        /**
         * Get bounds into our vectors given the new weight and calculate the number of boxes that would result.
         */
        unsigned
        calculate_new_bounds( weight new_weight, iterators_by_char & new_bounds ) {
            unsigned total_num_boxes = 0;
            BOOST_FOREACH( typename seqs_by_char::value_type seqs_for_char, seqs ) {
                new_bounds[ seqs_for_char.first ].clear();
                unsigned num_boxes_for_char = 1;
                BOOST_FOREACH( const sequence & seq, seqs_for_char.second ) {
                    iterator new_bound = boost::upper_bound(
                        seq.make_range(),
                        new_weight,
                        weight_less_than()
                    );
                    new_bounds[ seqs_for_char.first ].push_back( new_bound );
                    num_boxes_for_char *= seq.values->end() - new_bound;
                }
                total_num_boxes += num_boxes_for_char;
            }
            return total_num_boxes;
        }


        /**
         * Update the bounds with the new ones.
         */
        void
        update_with_new_bounds( weight new_weight, const iterators_by_char & new_bounds, unsigned num_boxes ) {
            const bool update_lower = num_boxes >= max_num_boxes;
            if( update_lower ) {
                BOOST_ASSERT( lower_weight != new_weight );
                lower_weight = new_weight;
            } else {
                BOOST_ASSERT( upper_weight != new_weight );
                upper_weight = new_weight;
            }
            BOOST_FOREACH( typename iterators_by_char::value_type bounds_for_char, new_bounds ) {
                BOOST_ASSERT( bounds_for_char.second.size() == seqs[ bounds_for_char.first ].size() );
                for( unsigned i = 0; bounds_for_char.second.size() != i; ++i ) {
                    if( update_lower ) {
                        seqs[ bounds_for_char.first ][ i ].lower = bounds_for_char.second[ i ];
                    } else {
                        seqs[ bounds_for_char.first ][ i ].upper = bounds_for_char.second[ i ];
                    }
                }
            }
        }


        /**
         * Recursively bring the upper and lower bounds closer whilst keeping both on opposite sides
         * of the threshold that results in the correct number of boxes.
         */
        void
        restrict_bounds() {

            while( true ) { // until the bounds are close enough

                check_invariants();

                //
                // How close are the bounds together?
                //
                const weight weight_range = upper_weight - lower_weight;
                BOOST_ASSERT( weight_range > 0. );
                if( weight_range < 1e-5 * ( max_weight - min_weight ) ) {
                    break; // we are close enough
                }

                //
                // The new weight to test is half-way in-between the upper and lower.
                //
                BOOST_ASSERT( lower_weight < upper_weight );
                const weight new_weight = lower_weight + weight_range / 2.;
                if( new_weight >= upper_weight || new_weight <= lower_weight ) { // possible rounding errors
                    break; // cannot go any further due to machine precision or no range
                }

                //
                // Update the bounds based on the new weight.
                //
                iterators_by_char new_bounds;
                const unsigned new_num_boxes = calculate_new_bounds( new_weight, new_bounds );
                update_with_new_bounds( new_weight, new_bounds, new_num_boxes );
            }

//                BOOST_FOREACH( typename seqs_by_char::value_type seqs_for_char, seqs ) {
//                    BOOST_FOREACH( const sequence & seq, seqs_for_char.second ) {
//                        std::cerr
//                            << "'" << int( seqs_for_char.first ) << "' : Lower has "
//                            << seq.values->end() - seq.lower
//                            << ", Upper has " << seq.values->end() - seq.upper << " values\n";
//                    }
//                }
        }


        /**
         * Remove values such that we do not over-step the limit on the number of boxes.
         */
        void
        reduce_values() {

            using namespace boost;
            using boost::lambda::_1;

            //
            // Calculate how many boxes for each character.
            //
            char_count_vec num_boxes_per_char;
            BOOST_FOREACH( typename seqs_by_char::value_type seqs_for_char, seqs ) {
                const unsigned num_boxes_lower = calculate_num_boxes( seqs_for_char.second, false );
                const unsigned num_boxes_upper = calculate_num_boxes( seqs_for_char.second, true  );
                num_boxes_per_char.push_back(
                    char_count(
                        num_boxes_lower, // num boxes allowed for this character
                        character( seqs_for_char.first ), // character
                        num_boxes_upper, // num at upper bound
                        num_boxes_lower  // num at lower bound
                    )
                );
            }
            // sort, largest first, i.e. those characters for which we can remove the most boxes between lower and upper
            sort( num_boxes_per_char | adaptors::reversed, counts_less_than() );

            //
            // Whilst we have too many boxes
            //
            unsigned total = std::numeric_limits< unsigned >::max();
            while( true ) {
                //
                // check we are still sorted
                //
//                std::cerr << num_boxes_per_char.size() << " : ";
//                copy( num_boxes_per_char, std::ostream_iterator< char_count >( std::cerr, " " ) );
//                std::cout << "\n";
                BOOST_ASSERT( is_sorted( num_boxes_per_char | adaptors::reversed, counts_less_than() ) );

                //
                // How many boxes do we have at the lower bound?
                //
                total = accumulate(
                    num_boxes_per_char | adaptors::transformed( get_count() ),
                    unsigned( 0 )
                );
                if( total > max_num_boxes ) {
                    unsigned to_remove = total - max_num_boxes;

                    //
                    // The counts are sorted in descending order of how many can be removed from them. We will lower just those counts that
                    // have the most boxes in-between lower and upper.
                    //
                    BOOST_ASSERT( ! num_boxes_per_char.empty() );
                    char_count_it i = num_boxes_per_char.begin();
                    const unsigned most_boxes = can_remove()( *i );
                    char_count_it next_highest = find_if(
                        num_boxes_per_char,
                        lambda::bind( can_remove(), _1 ) != most_boxes
                    );
                    unsigned counts_to_reduce = next_highest - num_boxes_per_char.begin();
                    //
                    // Make sure we don't take so many away that the front of our count vec
                    // has lower values than the next highest
                    //
                    if( num_boxes_per_char.end() != next_highest ) {
                        BOOST_ASSERT( i->get< 0 >() > next_highest->get< 0 >() );
                        to_remove = std::min(
                            ( can_remove()( *i ) - can_remove()( *next_highest ) ) * counts_to_reduce,
                            to_remove
                        );
                    }
                    for( ; i != next_highest; ++i, --counts_to_reduce ) {
                        const unsigned removing = to_remove / counts_to_reduce;
                        i->get< 0 >() -= removing;
                        to_remove -= removing;
                        // check we didn't remove too much from the counts
                        BOOST_ASSERT( i->get< 2 >() <= i->get< 0 >() );
                    }
                    BOOST_ASSERT( 0 == counts_to_reduce );

                    //
                    // Make sure front of our counts vector is sorted
                    //
                    BOOST_ASSERT(
                        is_sorted(
                            make_iterator_range( num_boxes_per_char.begin(), next_highest ) | adaptors::reversed,
                            counts_less_than()
                        )
                    );
//                    sort(
//                        make_iterator_range( num_boxes_per_char.begin(), next_highest ) | adaptors::reversed,
//                        counts_less_than()
//                    );

                } else {
                    BOOST_ASSERT( total == max_num_boxes );
                    break;
                }
            }

            //
            // So now we now how many boxes each character should provide to limit the total, we need
            // to reduce the number of values for each character.
            //
            BOOST_FOREACH( const char_count & count, num_boxes_per_char ) {
                using std::log;
                using std::exp;

                const character c = count.template get< 1 >();
                sequence_vec & sequences = seqs[ c ];
                const double a = log( count.template get< 0 >() ); // log desired num boxes
                const double ldot = log( count.template get< 3 >() ); // log num boxes at lower threshold
                const double udot = log( std::max( unsigned( sequences.size() ), count.template get< 2 >() ) ); // log num boxes at upper threshold - must keep at least one per sequence
                BOOST_ASSERT( ldot >= a    );
                BOOST_ASSERT( a    >  udot );
                if( ldot > a ) {

                    //
                    // for each sequence of values
                    //
                    BOOST_FOREACH( const sequence & seq, sequences ) {

                        //
                        // Shuffle the values in between the lower and upper bounds to
                        // try to remove any bias we might have when removing values
                        //
                        random_shuffle(
                            *seq.values
                            | adaptors::sliced(
                                seq.lower - seq.values->begin(),
                                seq.upper - seq.values->begin()
                            )
                        );

                        //
                        // Remove the desired proportion of values
                        //
                        if( seq.upper != seq.lower ) {
                            const double ls = std::log( seq.values->end() - seq.lower );
                            const double us = std::log( std::max( 1, int( seq.values->end() - seq.upper ) ) ); // we need at least one value in each sequence
                            const double astars = ( ls - us ) / ( ldot - udot ) * ( a - udot ) + us;
                            const unsigned num_to_keep = exp( astars );
                            BOOST_ASSERT( num_to_keep > 0 );
                            typename value_vec::iterator erase_up_to = seq.values->end() - num_to_keep;
                            BOOST_ASSERT( erase_up_to <= seq.upper );
                            BOOST_ASSERT( seq.lower <= erase_up_to );
                            seq.values->erase( seq.values->begin(), erase_up_to );
                        } else {
                            // do erase up to lower though as we definitely want to ignore every value worse than the lower bound
                            seq.values->erase( seq.values->begin(), seq.values->begin() + ( seq.lower - seq.values->begin() ) );
                        }
                    }

                    //
                    // Check the number of boxes is correct
                    //
                    BOOST_ASSERT(
                        boost::accumulate(
                            sequences
                                | adaptors::transformed( lambda::bind( &sequence::values, _1 ) )
                                | adaptors::indirected
                                | adaptors::transformed( get_size< value_vec >() ),
                            unsigned( 1 ),
                            std::multiplies< unsigned >()
                        ) <= count.template get< 0 >()
                    );
                }
            }
        }
    };


    /**
     * Calculate the boxes for all the values indexed by character.
     */
    template<
        typename BoxTraits,    /**< Our output type traits. */
        typename OutIt          /**< Output iterator to write boxes to. */
    >
    static
    void
    get_boxes_for(
        const map_by_char & values_by_char,
        OutIt output_it
    ) {
        BOOST_FOREACH( const typename map_by_char::value_type & values, values_by_char ) {
            get_boxes_for< BoxTraits >( *( values.second ), output_it );
        }
    }


    /**
     * Calculate the boxes for all the values for one particular character.
     */
    template<
        typename BoxTraits,    /**< Our output type traits. */
        typename OutIt          /**< Output iterator to write boxes to. */
    >
    static
    void
    get_boxes_for(
        const value_vec_vec & values,
        OutIt output_it
    ) {
        using namespace boost;

        typedef typename BoxTraits::weight box_weight;
        typedef typename BoxTraits::data data;
        typedef typename BoxTraits::coord_array point;

        //
        // Check the box weight type matches the values weight type.
        //
        typedef is_same< weight, box_weight > are_weights_same_type;
        BOOST_STATIC_ASSERT_MSG( are_weights_same_type::value, "Weight types should be equal." );

        //
        // Check that the dimensions of our boxes match the size of the values collection
        //
        const unsigned dimensions = size( values );
        if( point::static_size != dimensions )
        {
            throw std::invalid_argument(
                BIO_MAKE_STRING( "Wrong # of input sequences: " << unsigned( point::static_size ) << " != " << dimensions ) );
        }

        //
        // Iterate over every combination of values from each sequence
        //
        for( multi_iterator it( values ); ! it.at_end(); it.next() ) {

            //
            // calculate the box for this combination
            //
            point start;
            point end;
            weight w = weight( 0.0 );
            for( unsigned d = 0; dimensions != d; ++d ) { // for each dimension

                //
                // fill in the values for this dimension
                //
                const value & v = values[ d ][ it._ind[ d ] ];
                start[ d ] = ValueTraits::get_start( v );
                end[ d ] = ValueTraits::get_end( v );
                w += ValueTraits::get_weight( v );
            }
            w /= dimensions; // average the weight over the dimensions

            //
            // Output the box through the output iterator
            //
            *output_it =
                BoxTraits::make_box(
                    ValueTraits::get_data( values[ 0 ][ it._ind[ 0 ] ] ), // get the data from the first value...
                    w,
                    start,
                    end
                );
            ++output_it;
        }
    }
};




BIO_NS_END

#endif //BIO_MAX_CHAIN_BOXES_H_

