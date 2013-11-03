/**
@file

Copyright John Reid 2006

*/

#include "biopsy/defs.h"
#include "biopsy/binding_hits_max_chain.h"
#include "biopsy/analyse.h"
#include "biopsy/pssm.h"

#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/indirected.hpp>
#include <boost/range/algorithm/copy.hpp>

namespace biopsy {

binding_hit_traits::data
binding_hit_traits::get_data( const binding_hit & h )
{
    return boost::addressof( h );
}

const binding_hit_traits::character &
binding_hit_traits::get_char( const binding_hit & h )
{
    return h._binder_name;
}

binding_hit_traits::weight
binding_hit_traits::get_weight( const binding_hit & h )
{
    const weight odds_ratio = get_odds_ratio_from_p_binding( h._p_binding );
    const weight pre_prior_odds_ratio = odds_ratio / pssm_parameters::singleton().binding_background_odds_prior;
    return ( pre_prior_odds_ratio < 1.0 ) ? 0.0 : log( pre_prior_odds_ratio );
}

binding_hit_traits::coord
binding_hit_traits::get_start( const binding_hit & h )
{
    return h._location._position;
}

binding_hit_traits::coord
binding_hit_traits::get_end( const binding_hit & h )
{
    return h._location._position + h._location._length;
}

///**
//Estimate the number of boxes the maximal chain algorithm will use at the given threshold.
//*/
//unsigned
//num_boxes_for_threshold(
//    binding_hits_vec_ptr hit_array,
//    double threshold );
//
///**
//Calculate the best threshold to filter hits at before calculating maximal chain.
//*/
//double
//calculate_best_threshold(
//    binding_hits_vec_ptr hit_array,
//    unsigned num_boxes_limit );
//

unsigned
get_max_chain_max_num_sequences()
{
    return BIO_MAX_CHAIN_MAX_SEQUENCES;
}


unsigned
num_hits_at_threshold( const binding_hit::vec & hits, double threshold ) {
    return std::count_if(
        hits.begin(),
        hits.end(),
        binding_hit::p_binding_greater( threshold )
    );
}


unsigned
num_boxes_for_threshold(
    binding_hits_vec_ptr hit_array,
    double threshold )
{
    unsigned num_boxes = 1;
    BOOST_FOREACH( const binding_hit::vec_ptr & hits, *hit_array )
    {
        num_boxes *= num_hits_at_threshold( *hits, threshold );
    }

    return num_boxes;
}


typedef std::pair< double, unsigned > size_for_threshold;
typedef std::pair< size_for_threshold, size_for_threshold > threshold_pair;


threshold_pair
calculate_best_threshold(
    binding_hits_vec_ptr hit_array,
    unsigned num_boxes_limit )
{
    //std::cout << "Calculating best threshold..." << std::endl;

    size_for_threshold upper( 1.0, num_boxes_for_threshold( hit_array, 1.0 ) );
    size_for_threshold lower( 0.0, num_boxes_for_threshold( hit_array, 0.0 ) );

    // we only have to do something if the lower limit is too large...
    if( lower.second > num_boxes_limit ) {

        while( upper.first - lower.first > 1e-5 ) //whilst our search has not converged sufficiently
        {
            BOOST_ASSERT( upper.second <= num_boxes_limit );
            BOOST_ASSERT( lower.second > num_boxes_limit );

            // the mid-point of the current thresholds
            const double next_threshold = ( upper.first + lower.first ) / 2.0;
            // how many boxes for the new threshold?
            const size_for_threshold next( next_threshold, num_boxes_for_threshold( hit_array, next_threshold ) );

            // which threshold to update, the lower or the upper?
            if( next.second > num_boxes_limit ) //do we have too many boxes still?
            {
                lower = next; //yes - so update lower bound
            }
            else
            {
                upper = next; //no - so update upper bound
            }
        }

        BOOST_ASSERT( lower.second > num_boxes_limit );
    }
    BOOST_ASSERT( upper.second <= num_boxes_limit );

    return threshold_pair( lower, upper );
}


/** Calculate how many to remove given a count and a log-scale reduction. */
unsigned
how_many_to_remove( unsigned count, double log_to_remove_per_count ) {
    const double log_count = std::log( count );
    const unsigned to_remove = std::ceil( count - std::exp( log_count - log_to_remove_per_count ) );
    const unsigned result = std::min( count - 1, std::max( 0u, to_remove ) );
    return result;
}


/**
 * We know how many hits we want to remove between which thresholds.
 */
binding_hit::vec_ptr
remove_thresholded_hits(
    const binding_hit::vec & hits,
    double lower,
    double upper,
    unsigned to_remove
) {
    BOOST_ASSERT( to_remove ); // we should want to remove at least 1 to call this function.

    // how many hits in total do we have to remove from?
    const unsigned to_remove_from_total = num_hits_at_threshold( hits, lower ) - num_hits_at_threshold( hits, upper );
    unsigned could_remove_from = 0;
    unsigned removed = 0;

    binding_hit::p_binding_greater lower_pred( lower );
    binding_hit::p_binding_greater upper_pred( upper );

    binding_hit::vec_ptr results( new binding_hit::vec );
    BOOST_FOREACH( const binding_hit & hit, hits ) {
        if( upper_pred( hit ) ) { // if we pass the upper threshold we always keep the hit
            results->push_back( hit );
        } else if( lower_pred( hit ) ) { // if between the lower and upper thresholds we may/may not
            if( removed * to_remove_from_total > could_remove_from * to_remove ) {
                results->push_back( hit ); // we keep this hit
            } else {
                ++removed; // we removed this hit
            }
            ++could_remove_from;
        }
    }
    BOOST_ASSERT( could_remove_from == to_remove_from_total );
#ifndef NDEBUG
    const unsigned should_be_left =
#endif
    to_remove_from_total - to_remove;
#ifndef NDEBUG
    const unsigned are_left =
#endif
    num_hits_at_threshold( *results, lower ) - num_hits_at_threshold( *results, upper );
    BOOST_ASSERT( should_be_left == are_left );

    return results;
}


namespace { //anonymous

struct max_chain_builder
{
    binding_hit::vec_ptr _mc;
    max_chain_builder( binding_hit::vec_ptr mc ) : _mc( mc ) { }
    template< typename BoxPtr >
    void operator()( BoxPtr b )
    {
        _mc->push_back( *( b->_data ) );
    }
};

} //anonymous namespace




binding_hit::vec_ptr
analyse_max_chain(
    binding_hits_vec_ptr hit_array,
    unsigned num_boxes_limit
) {
    using namespace boost;

    //calculate the maximal chain if we don't have too many sequences
    binding_hit::vec_ptr mc;
    if( BIO_MAX_CHAIN_MAX_SEQUENCES >= hit_array->size() && hit_array->size() != 1 )
    {
        //calculate the maximal chain if possible
        mc.reset( new binding_hit::vec );
        if(
            ! BIO_NS::max_chain< binding_hit_traits >(
                *hit_array | adaptors::indirected,
                make_function_output_iterator( max_chain_builder( mc ) ),
                num_boxes_limit
            )
        ) {
            mc.reset();
        }
    }

    return mc;
}



} //namespace biopsy

