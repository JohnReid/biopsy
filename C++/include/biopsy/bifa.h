/**
 * Copyright John Reid 2010
 *
 * \file Implement BiFa algorithm.
 */


#ifndef BIOPSY_BIFA_H_
#define BIOPSY_BIFA_H_


#include <biopsy/defs.h>
#include <boost/range.hpp>
#include <boost/multi_array.hpp>
#include <algorithm>

namespace biopsy {
namespace bifa {




/// Convert a char to an int representation of a base
struct convert_char_base_to_int {

	/// Result type is int.
	typedef int result_type;

	/// Convert a char base to an int.
	int operator()( char b ) const {
		switch( b ) {
		case 'a': case 'A': return 0;
		case 'c': case 'C': return 1;
		case 'g': case 'G': return 2;
		case 't': case 'T': return 3;
		case 'n': case 'N': return 4;
		}
		throw std::logic_error( BIOPSY_MAKE_STRING( "Unknown base in sequence: \"" << b << "\"" ) );
	}
};



/// Defines PSSM traits.
template< typename PSSM >
struct PssmTraits {
	typedef typename boost::range_reference< PSSM >::type        column;                   ///< The type of a column of the PSSM.
	typedef typename boost::range_difference< PSSM >::type       size_t;                   ///< The size type.
    typedef typename PSSM::template array_view< 2 >::type       reverse_complement;       ///< The type of the reverse complement of the PSSM.
};



/// DNA alphabet
struct DnaAlphabet {
	enum {
		unknown_char = 4
	};
};





/// Score a PSSM on a word
template<
	typename Alphabet,
	typename PSSM,
	typename SeqIt
>
double
score_word(
	PSSM            pssm,
	SeqIt           word_begin
) {
	double result = 0.;
	BOOST_FOREACH( auto column, pssm ) {

		// is the base known?
		if( Alphabet::unknown_char == *word_begin ) {

			// unknown base so score as worst
			result += *std::min_element( boost::begin( column ), boost::end( column ) );

		} else { // known base

			// add score to result
			result += column[ *word_begin ];

		}
		++word_begin;
	}

	return result;
}




/// Score a PSSM on one strand of a sequence
template<
	typename Alphabet,
	typename PSSM,
	typename Seq,
	typename OutIt
>
void
score_one_strand(
	PSSM            pssm,
	Seq             seq,
	OutIt           out_it
) {
	// do nothing if sequence not long enough
	if( boost::size( seq ) < boost::size( pssm ) ) {
		return;
	}

	// work out where start of last word is
	typedef typename boost::range_iterator< Seq >::type seq_it;
	seq_it seq_begin = boost::begin( seq );
	seq_it seq_end = boost::end( seq ) - boost::size( pssm ) + 1;

	// score each word in the sequence
	while( seq_begin != seq_end ) {
		*out_it = score_word< Alphabet >( pssm, seq_begin );
		++out_it;
		++seq_begin;
	}
}


/// Return a view on the PSSM that is a reverse complement
template< typename PSSM >
typename PssmTraits< PSSM >::reverse_complement
pssm_reverse_complement( PSSM & pssm ) {
    using boost::indices;
    using boost::size;
    typedef boost::multi_array_types::index_range range;
    const int alphabet_size = ( 0 == boost::size( pssm ) ) ? 0 : size( pssm[ 0 ] );

//  // cannot reverse complement an empty PSSM
//  if( 0 == boost::size( pssm ) ) {
//      throw std::logic_error( "Cannot make reverse complement of an empty PSSM." );
//        return
//            pssm[
//                indices
//                    [ range( 0, 0 ) ]
//                    [ range( 0, 0 ) ]
//            ];
//  }

    return pssm[ indices[ range( size( pssm )  - 1, -1, -1 ) ]
                         [ range( alphabet_size - 1, -1, -1 ) ] ];
}


/// Stores the likelihoods of the words in the sequence
struct sequence_likelihoods {
	typedef boost::shared_ptr< sequence_likelihoods > ptr; ///< Shared pointer type.

	/// Get the log likelihood of the given word in the sequence
	virtual double get_word_log_likelihood( size_t position, size_t word_length ) const = 0;
};


/// Uniform sequence likelihoods
struct uniform_sequence_likelihoods : sequence_likelihoods {
	virtual double get_word_log_likelihood( size_t position, size_t word_length ) const {
		return std::log( .25 ) * word_length;
	}
};




/// Score a PSSM (in log likelihood form) on both strands of a sequence
template<
	typename Alphabet,
	typename PSSM,
	typename BgLikelihoods,
	typename Seq,
	typename OutFn
>
void
score_sequence(
	PSSM            pssm,
	BgLikelihoods   bg_likelihoods,
	Seq             seq,
	OutFn           out_fn
) {
	typedef PssmTraits< PSSM > pssm_traits;

	// do nothing if sequence not long enough
	const typename pssm_traits::size_t word_size = boost::size( pssm );
	if( typename pssm_traits::size_t( boost::size( seq ) ) < word_size ) {
		return;
	}

	// work out where start of last word is
	typedef typename boost::range_iterator< Seq >::type seq_it;
	seq_it seq_end = boost::end( seq ) - word_size + 1;

	// score each word on the sequence's positive strand
	{
		seq_it seq_begin = boost::begin( seq );
		for( size_t pos = 0; seq_begin != seq_end; ++pos, ++seq_begin ) {
			out_fn( score_word< Alphabet >( pssm, seq_begin ) - bg_likelihoods.get_word_log_likelihood( pos, word_size ), pos, true );
		}
	}

	// score each word on the sequence's negative strand
	{
		typename pssm_traits::reverse_complement pssm_rev_comp = pssm_reverse_complement( pssm );
		seq_it seq_begin = boost::begin( seq );
		for( size_t pos = 0; seq_begin != seq_end; ++pos, ++seq_begin ) {
			out_fn( score_word< Alphabet >( pssm_rev_comp, seq_begin ) - bg_likelihoods.get_word_log_likelihood( pos, word_size ), pos, false );
		}
	}
}



} // namespace bifa
} // namespace biopsy


#endif //BIOPSY_BIFA_H_
