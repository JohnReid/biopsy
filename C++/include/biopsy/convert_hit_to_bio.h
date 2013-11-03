/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_CONVERT_HIT_TO_BIO_H_
#define BIOPSY_CONVERT_HIT_TO_BIO_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"

#include <bio/run_match.h>

namespace biopsy {



struct transform_biopsy_hit_to_bio
	: std::unary_function< biopsy::binding_hit, BIO_NS::MatchResults >
{
	BIO_NS::MatchResults operator()( const biopsy::binding_hit & hit ) const
	{
		USING_BIO_NS;

		return
			MatchResults( 
				parse_table_link_accession_number( hit._binder_name ),
				Hit(
					BIO_NS::float_t( hit._p_binding ),
					hit._location._position,
					! hit._location._positive_strand ) );
	}
};


template< 
	typename HitRange,
	typename OutputIt
>
void
convert_biopsy_hit_to_bio(
	const HitRange & biopsy_hits,
	OutputIt output_it )
{
	BOOST_FOREACH( const biopsy::binding_hit & hit, biopsy_hits )
	{
		*output_it++ = transform_biopsy_hit_to_bio()( hit );
	}
}



} //namespace biopsy



#endif //BIOPSY_CONVERT_HIT_TO_BIO_H_
