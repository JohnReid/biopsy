/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/run_match.h"

#include <algorithm>

#include <gsl/gsl_errno.h>




BIO_NS_START


MatchResults::MatchResults(
	TableLink link,
	Hit result)
	: link(link)
	, result(result)
{
}

MatchResults::MatchResults()
{
}

void sort_by_position(match_result_vec_t & matches)
{
	std::sort(matches.begin(), matches.end(), MatchResultsPositionLessThan());
}


RaiseHitToPower::RaiseHitToPower(double power)
	: power(power)
{
}

void RaiseHitToPower::operator()(MatchResults & hit) const
{
	gsl_sf_result gsl_result;

	//ignore zeros
	if ( hit.result.score != 0.0 )
	{

		if ( GSL_SUCCESS != gsl_sf_log_e( double( hit.result.score ), &gsl_result ) )
		{
			throw BIO_MAKE_STRING( "Could not take logarithm of " << hit.result.score );
		}

		if ( GSL_SUCCESS != gsl_sf_exp_e( power * gsl_result.val, &gsl_result ) )
		{
			throw BIO_MAKE_STRING( "Could not take exponent of " << power * gsl_result.val );
		}

		hit.result.score = float_t( gsl_result.val );

	}
}


HitAboveThreshold::HitAboveThreshold(float_t threshold)
	: threshold(threshold)
{
}

bool HitAboveThreshold::operator()(const MatchResults & results) const
{
	return results.result.score > threshold;
}



BIO_NS_END

