/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_SEQUENCE_H_
#define BIOPSY_SEQUENCE_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "biopsy/defs.h"




namespace biopsy {



/**
Gets the complement of a nucleotide.
*/
struct nucleo_complement
	: std::unary_function< char, char >
{
	char operator()( char c ) const;
};



/**
Is the base unknown?
*/
struct is_unknown_nucleotide
{
	bool operator()( char c ) const;
};


/**
Is the base known, i.e. 'a', 'c', 'g' and 't'?
*/
struct is_known_nucleotide
{
	bool operator()( char c ) const;
};


/**
Is the sequence made of 'a', 'c', 'g' and 't's?
*/
struct is_known_sequence
{
	bool operator()( const std::string & seq ) const;
};


/**
Generates a random sequence. The seed must be non-zero.
*/
sequence
generate_random_sequence( unsigned length, unsigned seed );

/**
Returns an iterator that iterates over the complement of the original iterator.
*/
template< typename Iterator >
boost::transform_iterator< nucleo_complement, Iterator >
make_complement_iterator( Iterator it )
{
	return boost::make_transform_iterator( it, nucleo_complement() );
}





/**
Returns an iterator that iterates over the reverse complement of the original iterator.
*/
template< typename BidirectionalIterator >
boost::reverse_iterator< boost::transform_iterator< nucleo_complement, BidirectionalIterator > >
make_reverse_complement_iterator( BidirectionalIterator x )
{
	return boost::make_reverse_iterator( make_complement_iterator( x ) );
}



std::string
reverse_complement( const std::string & s );



} //namespace biopsy

#endif //BIOPSY_SEQUENCE_H_

