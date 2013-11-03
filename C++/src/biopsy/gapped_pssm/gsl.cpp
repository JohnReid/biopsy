/**
@file

Copyright John Reid 2006

*/

#include "biopsy/gsl.h"

#include <gsl/gsl_rng.h>


namespace biopsy {

namespace detail {

static const gsl_rng_type * T = 0;
static gsl_rng * r = 0;

extern "C" {

void
error_handler(
	const char * reason,
	const char * file,
	int line,
	int gsl_errno )
{
	throw std::logic_error( BIOPSY_MAKE_STRING( "GSL error: " << reason << " : " << file << " : " << gsl_errno ) );
}

} //extern "C"

} //namespace detail

gsl_rng *
get_gsl_rng()
{
	using namespace detail;

	if( 0 == T )
	{
		gsl_set_error_handler( error_handler );
		gsl_rng_env_setup();
		T = gsl_rng_default;
		if( ! T ) throw std::logic_error( "Could not initialise GSL random number generator" );
		r = gsl_rng_alloc (T);
		if( ! r ) throw std::logic_error( "Could not initialise GSL random number generator" );
	}
	return r;
}

} //namespace biopsy

