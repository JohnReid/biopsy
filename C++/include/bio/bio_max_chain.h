/**
@file

Copyright John Reid 2006

*/

#ifndef BIO_BINDING_MAX_CHAIN_H_
#define BIO_BINDING_MAX_CHAIN_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER


#include "bio/defs.h"
#include "bio/run_match.h"
#include "bio/max_chain.h"

#include <boost/function_output_iterator.hpp>
#include <boost/function.hpp>

BIO_NS_START


/** 
Describes some traits of binding hits for use in max chain algorithms. 
*/
struct binding_hit_traits
{
	typedef MatchResults hit;
	typedef const MatchResults * data;
	typedef TableLink character;
	typedef double weight;
	typedef int coord;

	static data get_data( const hit & h );
	static const character & get_char( const hit & h );
	static weight get_weight( const hit & h );
	static coord get_start( const hit & h );
	static coord get_end( const hit & h );
};




template< 
	unsigned d
>
struct binding_hit_max_chain_algorithm
{
	typedef typename max_chain_ns::MaxChainBox< 
		d, 
		binding_hit_traits::data, 
		binding_hit_traits::coord, 
		binding_hit_traits::weight 
	> box;


	typedef typename max_chain_ns::MaxChain<
		d, 
		box, 
		typename box::traits 
	> algorithm;



	template<
		typename HitSequenceRange,
		typename OutputIt
	>
	void
	operator()(
		const HitSequenceRange & hits,
		OutputIt output_it ) const
	{
	    typedef typename box::traits box_traits;
		typename box_traits::box_vec boxes;
		auto values = box_generator< binding_hit_traits >::template boxes_from_sequences< box_traits >(
			hits,
			std::back_inserter( boxes )
		);
		algorithm alg( boxes );

		BOOST_FOREACH( typename algorithm::box_ptr b, alg.maximal_chain ) {
			*output_it = b;
			++output_it;
		}
	}


	/**
	Print the box to the stream.
	*/
	static
	void
	print_box( typename algorithm::box_ptr b, std::ostream & os = std::cout )
	{
		std::copy( box::traits::get_start( *b )->begin(), box::traits::get_start( *b )->end(), std::ostream_iterator< int >( os, "," ) );
		os << "; ";
		std::copy( box::traits::get_end( *b )->begin(), box::traits::get_end( *b )->end(), std::ostream_iterator< int >( os, "," ) );
		os << "; " << box::traits::get_weight( *b ) << "; " << *box::traits::get_data( *b ) << std::endl;
	}
};




/**
A macro for each branch of the switch on the problem dimension.
*/
#define BIO_binding_hit_max_chain_case_instance( n ) \
	case n: \
	{ \
		typedef typename binding_hit_max_chain_algorithm< n >::algorithm::box_ptr arg; \
		boost::function1< void, arg > fn = boost::lambda::bind< void >( \
			&binding_hit_max_chain_algorithm< n >::print_box,  \
			boost::lambda::_1, \
			boost::ref( std::cout ) ); \
		binding_hit_max_chain_algorithm< n >()(  \
			hits,  \
			boost::make_function_output_iterator( fn ) ); \
	} \
	break;




/**
Calculates the maximal chain that is a subsequence of each sequence of hits in the HitSequenceRange.
*/
template< 
	typename HitSequenceRange
>
void
binding_hit_max_chain(
	const HitSequenceRange & hits )
{
	using namespace boost::lambda;
//	using boost::lambda::_1;

	if( 5 < boost::size( hits ) )
	{
		throw 
			std::logic_error( 
				BIO_MAKE_STRING( 
					"Max chain algorithm not implemented for " << boost::size( hits ) << " sequences. Limit is " << 5 ) );
	}

	switch( boost::size( hits ) )
	{
	case 0:
	{
		typedef typename binding_hit_max_chain_algorithm< 0 >::algorithm::box_ptr arg;
		boost::function1< void, arg > fn = boost::lambda::bind< void >(
			&binding_hit_max_chain_algorithm< 0 >::print_box, 
			boost::lambda::_1,
			boost::ref( std::cout ) );
		binding_hit_max_chain_algorithm< 0 >()(
			hits,
			boost::make_function_output_iterator( fn ) );
	}
	break;
	//BIO_binding_hit_max_chain_case_instance(  0 )
	BIO_binding_hit_max_chain_case_instance(  1 )
	BIO_binding_hit_max_chain_case_instance(  2 )
	BIO_binding_hit_max_chain_case_instance(  3 )
	BIO_binding_hit_max_chain_case_instance(  4 )
	BIO_binding_hit_max_chain_case_instance(  5 )
#if 0 //can increase compilation time significantly
	BIO_binding_hit_max_chain_case_instance(  6 )
	BIO_binding_hit_max_chain_case_instance(  7 )
	BIO_binding_hit_max_chain_case_instance(  8 )
	BIO_binding_hit_max_chain_case_instance(  9 )
	BIO_binding_hit_max_chain_case_instance( 10 )
	BIO_binding_hit_max_chain_case_instance( 11 )
#endif
	default:
		throw std::invalid_argument( "Cannot calculate the maximal chain for that many sequences." );
	}
}

#undef BIO_binding_hit_max_chain_case_instance


/**
Get a set of all the pssm tablelinks in biobase.
*/
template<
	typename OutputIt
>
void
get_all_pssm_links_in_biobase(
	const BiobasePssmFilter & filter,
	OutputIt output_it )
{
	BOOST_FOREACH( Matrix::map_t::value_type v, get_matrices( filter ) )
	{
		*output_it = v.first;
		++output_it;
	}
	BOOST_FOREACH( Site::map_t::value_type v, get_sites( filter ) )
	{
		*output_it = v.first;
		++output_it;
	}
}

/**
Returns estimate that the pssm binds the sequence. 
*/
template<
	typename OutputIt
>
float_t
score_pssm_over_sequence(
	TableLink link,
	seq_t::const_iterator seq_begin,
	seq_t::const_iterator seq_end,
	const MatchParams & params,
	OutputIt output_it )
{
	switch ( link.table_id )
	{
	case MATRIX_DATA:
		return
			score_pssm(
				*BiobaseDb::singleton().get_entry< MATRIX_DATA >( link ),
				seq_begin,
				seq_end,
				params,
				output_it );

	case SITE_DATA:
		return
			score_pssm(
				*BiobaseDb::singleton().get_entry< SITE_DATA >( link ),
				seq_begin,
				seq_end,
				params,
				output_it );

	default:
		throw std::logic_error( "Unknown pssm type" );
	}
}

template<
	typename OutputIt
>
void
get_pssms_in_results(
	const match_result_vec_t & hits,
	OutputIt output_it )
{
	BOOST_FOREACH( const MatchResults & hit, hits )
	{
		*output_it = hit.link;
		++output_it;
	}
}


/**
A vector of hit vectors, one for each sequence.
*/
typedef std::vector< match_result_vec_t > match_results_array;

/** 
Score sequences and store the match results in an output iterator. Maintains a set of PSSMs so that we only
score PSSMs that appear in all sequences.
OutputIt should take values convertible from match_result_vec_t.
*/
template<
	typename SequenceRange
>
void
score_phylo_sequences(
	const SequenceRange & sequences,
	const BiobasePssmFilter & filter,
	ScoreAlgorithm algorithm,
	float_t threshold,
	match_results_array & hit_array )
{
	const MatchParams params( threshold, algorithm );

	//get a set of all PSSMs in biobase
	typedef std::set< TableLink > table_link_set;
	table_link_set pssm_links;
	get_all_pssm_links_in_biobase( filter, std::inserter( pssm_links, pssm_links.begin() ) );

	//estimate the binding probabilities for each pssm in the phylo sequences
	typedef std::map< TableLink, float_t > binding_prob_map;
	binding_prob_map binding_prob_estimates;
	BOOST_FOREACH( const TableLink & pssm_link, pssm_links )
	{
		binding_prob_estimates[ pssm_link ] = float_t( 1.0 );
	}

	//for each sequence
	hit_array.clear();
	bool is_first_sequence = true;
	BOOST_FOREACH( const typename boost::range_value< SequenceRange >::type & sequence, sequences )
	{
		//push_back a hit vector for this sequence
		hit_array.push_back( match_result_vec_t() );
		match_result_vec_t & hits = *( hit_array.rbegin() );
		BOOST_FOREACH( const TableLink & pssm_link, pssm_links )
		{
			const float_t binding_p =
				score_pssm_over_sequence(
					pssm_link,
					sequence.begin(),
					sequence.end(),
					params,
					std::back_inserter( hits ) );

			//update the binding prob if it is not the first sequence.
			if( ! is_first_sequence )
			{
				binding_prob_estimates[ pssm_link ] *= binding_p;
			}
		}

		//rebuild set of pssms we are interested in
		pssm_links.clear();
		get_pssms_in_results( hits, std::inserter( pssm_links, pssm_links.begin() ) );

		//it won't be the first sequence next time around
		is_first_sequence = false;
	}

	//calculate the maximal chain
	binding_hit_max_chain( hit_array );

	//adjust the hits for phylogenetic conservation
	if( ! boost::empty( sequences ) )
	{
		BOOST_FOREACH( MatchResults & hit, hit_array[ 0 ] )
		{
			//adjust by estimate that binds in the phylo sequences and then take to
			//power(1/n) where n = # sequences
			hit.result.score =
				exp( 
					log( hit.result.score * binding_prob_estimates.find( hit.link )->second )
					/ ( boost::size( sequences ) ) );
		}
	}
}





BIO_NS_END

#endif //BIO_BINDING_MAX_CHAIN_H_

