#ifndef COMPRESSED_INT_ARRAY_H_
#define COMPRESSED_INT_ARRAY_H_


#include <vector>
#include <stdexcept>



/**
An array of integers, each with num_bits of storage.

Implemented as an array of larger ints (of type storage_int)
into which the client's ints are packed.
*/
template< typename storage_int, unsigned num_bits >
struct compressed_int_array
{
	typedef std::vector< storage_int > storage_array;

	storage_array storage;
	unsigned size;

	compressed_int_array( unsigned size = 0 );

	void resize( unsigned new_size );

	template< typename return_int >
	return_int get( unsigned idx ) const;

	template< typename value_int >
	void set( unsigned idx, value_int value );

	static unsigned storage_idx( unsigned idx );
	static unsigned in_storage_idx( unsigned idx );
};


template< typename storage_int, unsigned num_bits >
compressed_int_array< storage_int, num_bits >::compressed_int_array( unsigned s )
: size( 0 )
{
	if( ( 8 * sizeof(storage_int) ) % num_bits != 0 )
	{
		throw std::logic_error( "Must be a whole number of num_bits in storage_int" );
	}

	resize( s );
}

template< typename storage_int, unsigned num_bits >
void
compressed_int_array< storage_int, num_bits >::resize( unsigned new_size )
{
	size = new_size;
	storage.resize( size > 0 ? (1 + storage_idx( size - 1 )) : 0 );
}

template< typename storage_int, unsigned num_bits >
template< typename return_int >
return_int
compressed_int_array< storage_int, num_bits >::get( unsigned idx ) const
{
	if( idx >= size )
	{
		throw std::logic_error( "Not that many entries in array" );
	}

	const storage_int s = storage[ storage_idx( idx ) ];
	const unsigned in_idx = in_storage_idx( idx );
	return 
		(s << (8 * sizeof(storage_int) - (in_idx + 1) * num_bits)) 
		>> ((8 * sizeof( storage_int )) - num_bits);
}

template< typename storage_int, unsigned num_bits >
template< typename value_int >
void
compressed_int_array< storage_int, num_bits >::set( unsigned idx, value_int value )
{
	if( idx >= size )
	{
		throw std::logic_error( "Not that many entries in array" );
	}

	if( (value >> num_bits) != 0 )
	{
		throw std::logic_error( "Cannot store an integer that large" );
	}

	storage_int & s = storage[ storage_idx( idx ) ];
	const unsigned in_idx = in_storage_idx( idx );

	s ^= (value << (in_idx * num_bits));
}

template< typename storage_int, unsigned num_bits >
unsigned
compressed_int_array< storage_int, num_bits >::storage_idx( unsigned idx )
{
	return (idx * num_bits) / (8 * sizeof(storage_int));
}

template< typename storage_int, unsigned num_bits >
unsigned
compressed_int_array< storage_int, num_bits >::in_storage_idx( unsigned idx )
{
	return idx % ( (8 * sizeof(storage_int)) / num_bits );
}

#endif //COMPRESSED_INT_ARRAY_H_

