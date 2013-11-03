/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_LCS_H_
#define BIOPSY_LCS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "biopsy/defs.h"
#include "biopsy/binding_hits.h"

#include <bio/lcs.h>


namespace biopsy {



namespace detail {
struct binding_hit_char_extractor
{
	const std::string & operator()( const binding_hit & hit ) const
	{
		return hit._binder_name;
	}
};

struct binding_hit_start_extractor
{
	int operator()( const binding_hit & hit ) const
	{
		return hit._location._position;
	}
};

struct binding_hit_end_extractor
{
	int operator()( const binding_hit & hit ) const
	{
		return hit._location.get_end();
	}
};

struct binding_hit_score_extractor
{
	double operator()( const binding_hit & hit ) const
	{
		return hit._p_binding;
	}
};
} //namespace detail


typedef
	BIO_NS::LCS<
		binding_hit,
		std::string,
		detail::binding_hit_char_extractor,
		detail::binding_hit_start_extractor,
		detail::binding_hit_end_extractor,
		detail::binding_hit_score_extractor >
	lcs;

typedef boost::shared_ptr< lcs > lcs_ptr;

lcs_ptr lcs_create( const binding_hits_vec & hits_vec );
void lcs_calculate( lcs_ptr _lcs );
binding_hit::vec_ptr lcs_get( lcs_ptr _lcs );







} //namespace biopsy

#endif //BIOPSY_LCS_H_

