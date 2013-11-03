/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_BINDING_HITS_H_
#define BIOPSY_BINDING_HITS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"

namespace biopsy
{

	
class binding_hit_location
	: boost::equality_comparable< binding_hit_location >
	, boost::less_than_comparable< binding_hit_location >
{
public:
	int _position;
	bool _positive_strand;
	int _length;

public:
	binding_hit_location(
		int position = 0,
		int length = 0,
		bool positive_strand = true );

	bool operator==( const binding_hit_location & rhs ) const;
	bool operator<( const binding_hit_location & rhs ) const;

	int get_end() const;
};

std::ostream &
operator<<( std::ostream & os, const binding_hit_location & location );




class binding_hit
	: boost::equality_comparable< binding_hit >
	, boost::less_than_comparable< binding_hit_location >
{
public:
	std::string _binder_name;
	binding_hit_location _location;
	double _p_binding;

public:
	typedef std::vector< binding_hit > vec;
	typedef boost::shared_ptr< vec > vec_ptr;

public:
	binding_hit(
		const std::string & _binder_name = "",
		binding_hit_location _location = binding_hit_location(),
		double _p_binding = 0.0 );

	bool operator==( const binding_hit & rhs ) const;
	bool operator<( const binding_hit & rhs ) const;

	struct p_binding_greater
	{
		double threshold;
		p_binding_greater( double threshold );
		bool operator()( const binding_hit & hit ) const;
	};
};

typedef std::vector< binding_hit::vec_ptr > binding_hits_vec;
typedef boost::shared_ptr< binding_hits_vec > binding_hits_vec_ptr;

std::ostream &
operator<<( std::ostream & os, const binding_hit & hit );


/**
Returns a deduplicated vector of binder names that feature amongst the hits.
*/
string_vec_ptr
get_binder_names( binding_hit::vec_ptr hits );


/**
Sorts the hits by their positions.
*/
binding_hit::vec_ptr 
sort_hits_by_position( binding_hit::vec_ptr hits );



/**
An analysis of more than one sequence indexed by string.
*/
struct analysis
{
	typedef boost::shared_ptr< analysis > ptr;
	typedef std::string key;
	typedef std::map< key, binding_hit::vec_ptr > analysis_map;
	typedef std::vector< key > key_vec;
	typedef boost::shared_ptr< key_vec > key_vec_ptr;

	analysis_map _analyses;

	analysis();

	key_vec_ptr get_keys() const;
	binding_hit::vec_ptr get_hits_for( const key & k ) const;
	void set_hits_for( const key & k, binding_hit::vec_ptr hits );

	void serialise( const std::string & filename ) const;
	static analysis::ptr deserialise( const std::string & filename );
};



} //namespace biopsy

#endif //BIOPSY_BINDING_HITS_H_
