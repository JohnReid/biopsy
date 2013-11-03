/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_LOG_PROBS_H_
#define BIOPSY_LOG_PROBS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_machine.h>



#include <boost/test/floating_point_comparison.hpp>



namespace biopsy {


inline double safe_log( double x )
{
	return x <= 0.0
		? std::numeric_limits< double >::quiet_NaN()
		: gsl_sf_log( x );
}


template<
	typename ParameterRange
>
void
normalise_discrete_probabilities( 
	ParameterRange & parameters )
{
	using namespace boost;

	const double sum = std::accumulate( begin( parameters ), end( parameters ), 0.0 );
	BOOST_FOREACH( typename range_value< ParameterRange >::type & p, parameters )
	{
		p /= sum;
	}
}

template< 
	typename log_vector, 
	typename out_vector 
>
void
probabilities_from_logs_of_choose( 
	const log_vector & logs,
	out_vector & out )
{
	out.resize( logs.size() );

	// find the largest - as we can get underflow problems
	const double largest = *( std::max_element( logs.begin(), logs.end() ) );
	for( unsigned i = 0; logs.size() != i; ++i )
	{
		//scale by the largest
		const double v = logs[i] - largest;
		out[i] = (v <= GSL_LOG_DBL_MIN) ? 0.0 : gsl_sf_exp( v );
		BOOST_ASSERT( ! _isnan( out[i] ) );
	}

	normalise_discrete_probabilities( out );
}

template< 
	typename ValueRange, 
	typename ParameterRange
>
double
dirichlet_log_pdf( const ValueRange & values, const ParameterRange & parameters )
{
	using namespace boost;

	const double param_sum = std::accumulate( begin( parameters ), end( parameters ), 0.0 );

	double result = gsl_sf_lngamma( param_sum );

	typedef boost::tuple< double, double > param_value_tuple;
	BOOST_FOREACH(
		const param_value_tuple & t, 
		make_iterator_range(
			make_zip_iterator( make_tuple( begin( parameters ), begin( values ) ) ),
			make_zip_iterator( make_tuple( end( parameters ), end( values ) ) ) ) )
	{
		const double p = t.get<0>();
		const double v = t.get<1>();

		result -= gsl_sf_lngamma( p );
		result += ( p - 1.0 ) * gsl_sf_log( v );
	}

	return result;
}

template< 
	typename ParameterRange, 
	typename OutputRange
>
void
draw_from_dirichlet( const ParameterRange & alpha, OutputRange & output )
{
	using namespace boost;

	BOOST_ASSERT( size( alpha ) == size( output ) );

	gsl_ran_dirichlet(
		get_gsl_rng(),
		size( alpha ),
		&*begin( alpha ),
		&*begin( output ) );
}


template< 
	typename ParameterRange
>
unsigned
draw_from_multi( const ParameterRange & alpha )
{
	using namespace boost;
	using namespace boost::test_tools;

	BOOST_ASSERT( 
		close_at_tolerance< double >( percent_tolerance_t< double >( 1.0 ) )( 
			1.0,  
			std::accumulate( begin( alpha ), end( alpha ), 0.0 ) ) );

	double p = gsl_rng_uniform( get_gsl_rng() );

	unsigned result = 0;
	BOOST_FOREACH( const double a, alpha )
	{
		if( p < a )
		{
			break; 
		}
		p -= a;
		++result;
	}

	BOOST_ASSERT( result < size( alpha ) );

	return result;
}



} //namespace biopsy

#endif //BIOPSY_LOG_PROBS_H_
