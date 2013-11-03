/* Copyright John Reid 2007
*/

#include "bio-pch.h"


/**
@file

Copyright John Reid 2006

*/

#include "bio/bio_max_chain.h"
#include "bio/biobase_match.h"

using namespace boost;
using namespace std;

BIO_NS_START





binding_hit_traits::data
binding_hit_traits::get_data( const hit & h ) 
{
	return boost::addressof( h ); 
}

const binding_hit_traits::character & binding_hit_traits::get_char( const hit & h ) 
{ 
	return h.link; 
}

binding_hit_traits::weight 
binding_hit_traits::get_weight( const hit & h ) 
{
	return h.result.score; 
}

binding_hit_traits::coord 
binding_hit_traits::get_start( const hit & h ) 
{ 
	return h.result.position; 
}

binding_hit_traits::coord 
binding_hit_traits::get_end( const hit & h ) 
{ 
	return h.result.position + make_pssm( h.link ).size(); 
}




BIO_NS_END



