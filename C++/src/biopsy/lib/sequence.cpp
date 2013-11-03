/**
@file

Copyright John Reid 2006

*/

#include "biopsy/sequence.h"

using namespace boost;
using namespace std;

namespace biopsy {




char
nucleo_complement::operator()( char c ) const
{
	switch( c )
	{
	case 'a': return 't';
	case 'A': return 'T';
	case 'c': return 'g';
	case 'C': return 'G';
	case 'g': return 'c';
	case 'G': return 'C';
	case 't': return 'a';
	case 'T': return 'A';
	}
	throw std::logic_error( BIOPSY_MAKE_STRING( "Unknown nucleotide: '" << c << "'" ) );
}



std::string
reverse_complement( const std::string & s )
{
	std::string result;
	std::copy(
		make_reverse_complement_iterator( s.end() ),
		make_reverse_complement_iterator( s.begin() ),
		std::back_inserter( result ) );
	return result;
}


void
append_sequence( sequence_vec & v, const sequence & seq )
{
	v.push_back( seq );
}
void
append_sequences( sequence_vec & v, const sequence_vec & seqs )
{
	v.insert(v.end(),seqs.begin(),seqs.end() );
}


bool is_known_sequence::operator()( const std::string & seq ) const
{
	BOOST_FOREACH( char c, seq )
	{
		if( ! is_known_nucleotide()( c ) )
		{
			return false;
		}
	}
	return true;
}


bool is_known_nucleotide::operator()( char c ) const
{
	switch( c )
	{
	case 'a':
	case 'A':
	case 'c':
	case 'C':
	case 'g':
	case 'G':
	case 't':
	case 'T':
		return true;
	default:
		return false;
	}
}

bool is_unknown_nucleotide::operator()( char c ) const
{
	return c == 'n' || c == 'N';
}


namespace detail {

template< typename RNG >
struct
random_nucleotide
{
	RNG _rng;
	boost::uniform_int<> _four;
	boost::variate_generator< RNG, boost::uniform_int<> > _gen;

	random_nucleotide( RNG rng )
		: _rng( rng )
		, _four( 3 )
		, _gen( _rng, _four )
	{ }

	char operator()()
	{
		switch( _gen() )
		{
		case 0: return 'a';
		case 1: return 'c';
		case 2: return 'g';
		case 3: return 't';
		default: 
			BOOST_ASSERT( false );
			return ' ';
		}
	}
};

} //namespace detail


sequence
generate_random_sequence( unsigned length, unsigned s )
{
	typedef boost::mt19937 rng;
	rng my_rng;
	my_rng.seed( rng::result_type( s ) );
	detail::random_nucleotide< rng & > rand_nucleo( my_rng );

	sequence result( length, ' ' );
	std::generate( result.begin(), result.end(), rand_nucleo );
	return result;
}



} //namespace biopsy
