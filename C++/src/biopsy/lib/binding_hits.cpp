/**
@file

Copyright John Reid 2006

*/

#include "biopsy/defs.h"
#include "biopsy/binding_hits.h"

#include <bio/serialisable.h>


namespace boost {
namespace serialization {

template< typename Archive >
void serialize(
	Archive & archiver, 
	biopsy::binding_hit_location & archivee, 
	const unsigned int version )
{
    archiver & archivee._position;
    archiver & archivee._positive_strand;
    archiver & archivee._length;
}

template< typename Archive >
void serialize(
	Archive & archiver, 
	biopsy::binding_hit & archivee, 
	const unsigned int version )
{
    archiver & archivee._binder_name;
    archiver & archivee._location;
    archiver & archivee._p_binding;
}

template< typename Archive >
void serialize(
	Archive & archiver, 
	biopsy::analysis & archivee, 
	const unsigned int version )
{
    archiver & archivee._analyses;
}

} //namespace boost
} //namespace serialization



namespace biopsy
{

binding_hit_location::binding_hit_location(
	int position,
	int length,
	bool positive_strand )
	: _position( position )
	, _positive_strand( positive_strand )
	, _length( length )
{
}


bool
binding_hit_location::operator==( const binding_hit_location & rhs ) const
{
	return
		this == boost::addressof( rhs )
		|| (
			_position == rhs._position
			&&
			_positive_strand == rhs._positive_strand
			);
}

bool
binding_hit_location::operator<( const binding_hit_location & rhs ) const
{
	if( _position < rhs._position ) return true;
	if( rhs._position < _position ) return false;
	if( _positive_strand < rhs._positive_strand ) return true;
	return false;
}

int
binding_hit_location::get_end() const
{
	return _position + _length;
}


binding_hit::binding_hit(
	const std::string & binder_name,
	binding_hit_location location,
	double p_binding )
	: _binder_name( binder_name )
	, _location( location )
	, _p_binding( p_binding )
{
}


bool
binding_hit::operator==( const binding_hit & rhs ) const
{
	return
		this == boost::addressof( rhs )
		|| (
			_binder_name == rhs._binder_name
			&&
			_location == rhs._location
			);
}

bool
binding_hit::operator<( const binding_hit & rhs ) const
{
	if( _location < rhs._location ) return true;
	if( rhs._location < _location ) return false;
	if( _binder_name < rhs._binder_name ) return true;
	return false;
}

binding_hit::p_binding_greater::p_binding_greater( double threshold )
: threshold( threshold )
{
}

bool
binding_hit::p_binding_greater::operator()( const binding_hit & hit ) const
{
	return hit._p_binding >= threshold;
}

std::ostream &
operator<<( std::ostream & os, const binding_hit_location & location )
{
	os << location._position << ";" << location._length << ";" << ( location._positive_strand ? "+" : "-" );
	return os;
}


std::ostream &
operator<<( std::ostream & os, const binding_hit & hit )
{
	os << hit._binder_name << ";" << hit._p_binding << ";" << hit._location;
	return os;
}



string_vec_ptr
get_binder_names( binding_hit::vec_ptr hits )
{
	std::set< std::string > names;
	BOOST_FOREACH( const binding_hit & hit, *hits )
	{
		names.insert( hit._binder_name );
	}

	string_vec_ptr result( new string_vec );
	std::copy( names.begin(), names.end(), std::back_inserter( *result ) );

	return result;
}

struct hits_position_less_than
{
	bool operator()( const binding_hit & lhs, const binding_hit & rhs ) const
	{
		return lhs._location._position < rhs._location._position;
	}
};

binding_hit::vec_ptr 
sort_hits_by_position( binding_hit::vec_ptr hits )
{
	using namespace boost::lambda;
//	using boost::lambda::_1;

	std::sort(
		hits->begin(),
		hits->end(),
		hits_position_less_than() );

	return hits;
}


analysis::analysis()
{
}


analysis::key_vec_ptr
analysis::get_keys() const
{
	key_vec_ptr result( new key_vec );
	BOOST_FOREACH( analysis_map::value_type v, _analyses )
	{
		result->push_back( v.first );
	}
	return result;
}


binding_hit::vec_ptr
analysis::get_hits_for( const key & k ) const
{
	analysis_map::const_iterator i = _analyses.find( k );
	if( _analyses.end() == i )
	{
		throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find analysis for: " << k ) );
	}
	return i->second;
}

void
analysis::set_hits_for( const key & k, binding_hit::vec_ptr hits ) 
{
	_analyses[ k ] = hits;
}


void
analysis::serialise( const std::string & filename ) const
{
	BIO_NS::serialise< true >(
		*this,
		boost::filesystem::path( filename ) );
}

analysis::ptr
analysis::deserialise( const std::string & filename )
{
	return
		BIO_NS::deserialise< true, analysis >(
			boost::filesystem::path( filename ) );

}


} //namespace biopsy

