/**
@file

Copyright John Reid 2006

*/

#include <biopsy/gapped_pssm_hmm.h>
#include <biopsy/gsl.h>
#include <biopsy/log_probs.h>

#include <vector>
#include <sstream>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>

#include <limits>

#include <math.h>

#ifdef _WIN32
# define MY_ISNAN(x) _isnan(x)
#else
# define MY_ISNAN(x) isnan(x)
#endif

#define VERIFY( x ) if( ! ( x ) ) { throw std::logic_error( "Verify check failed" ); }
//#define VERIFY( x ) BOOST_ASSERT( x )

namespace biopsy {
namespace gapped_pssm {
namespace hmm {

void
calculate_all_transition_parameters(
	unsigned K,
	transition_parameter_array & parameter_array );

observed_data::observed_data(
	unsigned K,
	const observed_sequences & X,
	double_vector psi,
	double_vector theta,
	double_vector phi,
	double_vector upsilon )
	: X( X )
	, K( K )
	, Psi( psi )
	, Theta( theta )
{
	if( psi.size() != 4 ) throw std::logic_error( "psi must be length 4" );
	if( theta.size() != 4 ) throw std::logic_error( "theta must be length 4" );
	if( phi.size() != 2 ) throw std::logic_error( "phi must be length 2" );
	if( upsilon.size() != 2 ) throw std::logic_error( "upsilon must be length 2" );

	Phi[0] = phi[0];
	Phi[1] = phi[1];
	Upsilon[0] = upsilon[0];
	Upsilon[1] = upsilon[1];
}

hidden_data::hidden_data( const observed_data & data )
	: e( boost::extents[ data.E() ][ 4 ] )
	, t( data.K + 1 )
{
	calculate_all_transition_parameters( data.K, trans_params );

	for( unsigned e_idx = 0; data.E() != e_idx; ++e_idx ) 
	{
		draw_from_dirichlet( (0 == e_idx) ? data.Theta : data.Psi, e[ e_idx ] );
	}

	for( unsigned k_idx = 0; data.K + 1 != k_idx; ++k_idx )
	{
		t[ k_idx ] = gsl_ran_beta( 
			get_gsl_rng(),
			(0 == k_idx) ? data.Upsilon[0] : data.Phi[0],
			(0 == k_idx) ? data.Upsilon[1] : data.Phi[1] );
	}
}

void 
hidden_data::draw_sequence( unsigned I, unsigned_vector & r, sequence & seq ) const
{
	const state_map & map = get_state_map( t.size() - 1 ); //map for size K

	draw_r( I, r );

	seq.resize( I );
	for( unsigned i = 0; I != i; ++i )
	{
		const unsigned s = r[i];
		const unsigned k = map.k( s );
		const bool c = map.c( s );

		const unsigned drawn_x = draw_from_multi( e[ k ] );
		seq[ i ] = c ? 3 - drawn_x : drawn_x;
	}
}

void 
hidden_data::draw_r( unsigned I, unsigned_vector & r ) const
{
	r.resize( I, trans_params.size() );
	if( r.size() )
	{
		r[0] = 0; //first base is always background
		for( unsigned i = 1; I != i; ++i ) //for each base
		{
			double sample = gsl_rng_uniform( get_gsl_rng() );
			for( unsigned successor = 0; trans_params.size() != successor; ++successor ) //for each possible successor
			{
				//the parameters
				const transition_parameters & params = trans_params[ r[i-1] ][ successor ];

				if( params.get<0>() ) //if there is a transition
				{
					//the chance of this transition
					const double p_this_successor = params.get<1>() + params.get<2>() * t[ params.get<3>() ];

					if( sample < p_this_successor )
					{
						r[i] = successor;
						break;
					}
					else
					{
						sample -= p_this_successor;
					}
				}
			}
			VERIFY( r[i] != trans_params.size() );
		}
	}
}


variational_distribution::variational_distribution( const observed_data & data )
	: rho( data.X.N() )
	, eta( boost::extents[ data.E() ][ 4 ] )
	, tau( boost::extents[ data.K ][ 2 ] )
{
	for( unsigned n = 0; data.X.N() != n; ++n )
	{
		rho[n].resize( data.X.I(n), double_vector( data.S(), 1.0 / double( data.S() ) ) );
	}

	for( unsigned e = 0; data.E() != e; ++e ) 
	{
		draw_from_dirichlet( (0 == e) ? data.Theta : data.Psi, eta[ e ] );
	}

	for( unsigned k = 0; data.K != k; ++k ) 
	{
		tau[ k ][ 0 ] = ( k == 0 ) ? data.Upsilon[ 0 ] : data.Phi[ 0 ];
		tau[ k ][ 1 ] = ( k == 0 ) ? data.Upsilon[ 1 ] : data.Phi[ 1 ];
	}
}

void variational_distribution::r_mode( unsigned_vector_vec & r ) const
{
	r.resize( rho.size() );
	for( unsigned n = 0; r.size() != n; ++n )
	{
		r[n].resize( rho[n].size() );
		for( unsigned i = 0; r[n].size() != i; ++i )
		{
			r[n][i] = std::max_element( rho[n][i].begin(), rho[n][i].end() ) - rho[n][i].begin();
		}
	}
}

/**
Defines connectivity between states. I.e. given a state will call a output iterator once for each state that can preceed it.
*/
template< typename OutIt >
void generate_predecessor_states(
	const unsigned K,
	const unsigned s,
	OutIt out_it )
{
	const state_map & map = get_state_map( K );

	const bool g = map.g( s );
	const bool c = map.c( s );
	const unsigned k = map.k( s );

	//work out which states were possible predecessors
	if( map.b( s ) ) //is it background?
	{
		//background was possible predecessor
		*out_it++ = map.s( false, false, 0 );
		
		//non-rev-comp non-gap base K possible
		*out_it++ = map.s( false, false, K );
		
		//rev-comp non-gap base 1 possible
		*out_it++ = map.s( false, true, 1 );
	}
	else if( g ) //is it a gap base?
	{
		if( c ) //is it rev-comp?
		{
			//rev-comp non-gap base only possibility
			*out_it++ = map.s( false, true, k+1 );
		}
		else //not rev-comp
		{
			//non-rev-comp non-gap base only possibility
			*out_it++ = map.s( false, false, k );
		}
	}
	else //not a gap base
	{
		if( c ) //is it rev-comp?
		{
			//is it the first base?
			if( k == K )
			{
				//background only possible predecessor
				*out_it++ = map.s( false, false, 0 );
			}
			else //not the first base
			{
				//rev-comp gap base was possible predecessor
				*out_it++ = map.s( true, true, k );

				//rev-comp non-gap base was possible predecessor
				*out_it++ = map.s( false, true, k+1 );
			}
		}
		else //not rev-comp
		{
			//is it the first base?
			if( k == 1 )
			{
				//background only possible predecessor
				*out_it++ = map.s( false, false, 0 );
			}
			else
			{
				//non-rev-comp gap base was possible predecessor
				*out_it++ = map.s( true, false, k-1 );

				//non-rev-comp non-gap base was possible predecessor
				*out_it++ = map.s( false, false, k-1 );
			}
		}
	}
}

void calculate_predecessor_states( unsigned K, unsigned_vector_vec & predecessor_states )
{
	const unsigned S = get_state_map( K ).num_states();
	predecessor_states.resize( S );
	for( unsigned s = 0; S != s; ++s )
	{
		predecessor_states[ s ].clear();
		generate_predecessor_states( K, s, std::inserter( predecessor_states[ s ], predecessor_states[ s ].begin() ) );
	}
}



/** Build an array of transition probability parameters. */
void
calculate_all_transition_parameters(
	unsigned K,
	transition_parameter_array & parameter_array )
{
	using boost::make_tuple;

	const state_map & map = get_state_map( K );

	//make sure correct size
	parameter_array.resize( boost::extents[ map.num_states() ][ map.num_states() ] );

	//fill array with 0.0
	for( unsigned pre_s = 0; parameter_array.size() != pre_s; ++pre_s )
		for( unsigned s = 0; parameter_array[0].size() != s; ++s )
			parameter_array[ pre_s ][ s ] = boost::make_tuple( false, 0.0, 0.0, 0 );

	//
	//for each possible transition prob set value.
	//

	// background -> background
	parameter_array[ map.s( false, false, 0 ) ][ map.s( false, false, 0 ) ] = make_tuple( true, 1.0, -1.0, 0 );

	// background -> first of site (non-rev-comp and rev-comp)
	parameter_array[ map.s( false, false, 0 ) ][ map.s( false, false, 1 ) ] 
		= parameter_array[ map.s( false, false, 0 ) ][ map.s( false, true, K ) ] 
		= make_tuple( true, 0.0, 0.5, 0 );

	for( unsigned k = 1; K != k; ++k ) //we have K-1 transitions inside the pssm
	{
		// non-gap -> gap (non-rev-comp and rev-comp)
		parameter_array[ map.s( false, false, k ) ][ map.s( true, false, k ) ] 
			= parameter_array[ map.s( false, true, k+1 ) ][ map.s( true, true, k ) ] 
			= make_tuple( true, 0.0, 1.0, k );

		// non-gap -> non-gap (non-rev-comp and rev-comp)
		parameter_array[ map.s( false, false, k ) ][ map.s( false, false, k+1 ) ] 
			= parameter_array[ map.s( false, true, k+1 ) ][ map.s( false, true, k ) ] 
			= make_tuple( true, 1.0, -1.0, k );

		// gap -> non-gap (non-rev-comp and rev-comp)
		parameter_array[ map.s( true, false, k ) ][ map.s( false, false, k+1 ) ] 
			= parameter_array[ map.s( true, true, k ) ][ map.s( false, true, k ) ] 
			= make_tuple( true, 1.0, 0.0, k );
	}

	//last of site (non-rev-comp and rev-comp) -> background
	parameter_array[ map.s( false, false, K ) ][ map.s( false, false, 0 ) ]
		= parameter_array[ map.s( false, true, 1 ) ][ map.s( false, false, 0 ) ]
		= make_tuple( true, 1.0, 0.0, 0 );

}



void calculate_p_r_given_predecessor( 
	const transition_parameter_array & trans_params,
	const variational_distribution & var_dist,
	double_array & p_r_given_predecessor )
{
	//make sure correct sizes
	p_r_given_predecessor.resize( boost::extents[ trans_params.shape()[0] ][ trans_params.shape()[1] ] );

	//calculate array
	for( unsigned pre_s = 0; trans_params.shape()[0] != pre_s; ++pre_s )
	{
		for( unsigned s = 0; trans_params.shape()[1] != s; ++s )
		{
			const transition_parameters & params = trans_params[ pre_s ][ s ];

			p_r_given_predecessor[ pre_s ][ s ] =
				params.get<0>()
					? params.get<1>() + params.get<2>() * var_dist.expected_t( params.get<3>() )
					: 0.0;
		}
	}
}

void calculate_logs(
	const double_array & p,
	double_array & log_p )
{
	const double_array::size_type * shape = p.shape();
	log_p.resize( boost::extents[ shape[ 0 ] ][ shape[ 1 ] ] );
	for( unsigned i1 = 0; shape[0] != i1; ++i1 )
	{
		for( unsigned i2 = 0; shape[1] != i2; ++i2 )
		{
			log_p[ i1 ][ i2 ] = safe_log( p[ i1 ][ i2 ] );
		}
	}
}

void check_p_r_given_predecessor( const double_array & p_r_given_predecessor )
{
	for( unsigned pre_s = 0; p_r_given_predecessor.shape()[0] != pre_s; ++pre_s )
	{
		double sum = 0.0;
		for( unsigned s = 0; p_r_given_predecessor.shape()[1] != s; ++s )
		{
			sum += p_r_given_predecessor[ pre_s ][ s ];
		}
		VERIFY( boost::test_tools::check_is_close( 1.0, sum, boost::test_tools::percent_tolerance( 0.01f ) ) );
	}
}

void calculate_p_predecessor_given_r( 
	const double_array & p_r_given_predecessor,
	double_array & p_predecessor_given_r )
{
	const double_array::size_type * shape = p_r_given_predecessor.shape();
	p_predecessor_given_r.resize( boost::extents[ shape[ 0 ] ][ shape[ 1 ] ] );
	for( unsigned s = 0; shape[1] != s; ++s )
	{
		double sum = 0.0;
		for( unsigned pre_s = 0; shape[0] != pre_s; ++pre_s )
		{
			sum += p_r_given_predecessor[ pre_s ][ s ];
		}
		for( unsigned pre_s = 0; shape[0] != pre_s; ++pre_s )
		{
			p_predecessor_given_r[ pre_s ][ s ] = p_r_given_predecessor[ pre_s ][ s ] / sum;
		}
	}
}

void
calculate_log_p_x_given_r( unsigned K, unsigned S, const variational_distribution & var_dist, double_array & log_p_x_given_r )
{
	log_p_x_given_r.resize( boost::extents[ S ][ 4 ] );
	const state_map & map = get_state_map( K );
	for( unsigned s = 0; S != s; ++s )
	{
		const unsigned m = map.m( s );
		const bool c = map.c( s );
		for( unsigned x = 0; 4 != x; ++x )
		{
			log_p_x_given_r[ s ][ x ] = gsl_sf_log( var_dist.eta[ m ][ c ? 3 - x : x ] );
		}
	}
}

model::model( observed_data::ptr data )
: data( data )
, var_dist( new variational_distribution( *data ) )
{
	calculate_predecessor_states( data->K, predecessor_states );
	calculate_all_transition_parameters( data->K, trans_params );
}


double
model::update(
	bool update_rho, 
	bool update_eta, 
	bool update_tau )
{
	double LL = 0.0;

	const unsigned K = data->K;
	const unsigned S = data->S();
	const unsigned N = data->X.N();
	const unsigned E = data->E();
	const state_map & map = get_state_map( K );

	double_array p_r_given_predecessor;
	calculate_p_r_given_predecessor( trans_params, *var_dist, p_r_given_predecessor );

	double_array log_p_r_given_predecessor;
	calculate_logs( p_r_given_predecessor, log_p_r_given_predecessor );

	double_array p_predecessor_given_r;
	calculate_p_predecessor_given_r( p_r_given_predecessor, p_predecessor_given_r );

	double_array log_p_x_given_r;
	calculate_log_p_x_given_r( K, S, *var_dist, log_p_x_given_r );

	for( unsigned e = 0; E != e; ++e ) //for each base in the pssm and the background base
	{
		LL += dirichlet_log_pdf( var_dist->eta[ e ], (0 == e) ? data->Theta : data->Psi );
		VERIFY( ! MY_ISNAN( LL ) );

		if( update_eta )
		{
			for( unsigned x = 0; 4 != x; ++x )
			{
				var_dist->eta[ e ][ x ] = (0 == e) ? data->Theta[ x ] : data->Psi[ x ];
			}
		}
	}

	for( unsigned k = 0; K != k; ++k ) //for each state
	{
		LL += dirichlet_log_pdf( var_dist->tau[ k ], (0 == k) ? data->Upsilon : data->Phi );
		VERIFY( ! MY_ISNAN( LL ) );

		if( update_tau )
		{
			var_dist->tau[ k ][ 0 ] = (0 == k) ? data->Upsilon[ 0 ] : data->Phi[ 0 ];
			var_dist->tau[ k ][ 1 ] = (0 == k) ? data->Upsilon[ 1 ] : data->Phi[ 1 ];
		}
	}

	for( unsigned n = 0; N != n; ++n ) //for each sequence
	{
		const unsigned I = data->X.I( n );

		double_vector_vec log_rho( I ); //IxS
		if( update_rho )
		{
			for( unsigned i = 0; I != i; ++i ) //for each base
			{
				log_rho[ i ].resize( S, 0.0 );
			}
		}

		for( unsigned i = 0; I != i; ++i ) //for each base
		{
			const unsigned x = data->X._sequences[ n ][ i ];
			if( 4 == x ) continue; //ignore 'n's

			for( unsigned s = 0; S != s; ++s ) //for each state
			{
				const double p_s = var_dist->rho[ n ][ i ][ s ];
				const bool g_s = map.g( s );
				const bool c_s = map.c( s );
				const unsigned k_s = map.k( s );
				const unsigned m_s = map.m( s );
				const double log_p_x = log_p_x_given_r[ s ][ x ];

				if( update_rho ) 
				{
					log_rho[ i ][ s ] += log_p_x;
				}

				LL += log_p_x * p_s;
				VERIFY( ! MY_ISNAN( LL ) );

				if( update_eta ) 
				{
					var_dist->eta[ m_s ][ c_s ? 3 - x : x ] += p_s;
				}

				if( 0 != i ) //if we're not at the beginning of the sequence (because we can't look at transition to first base)
				{
					/**
					double_vector p_s_given_pre( S, 0.0 );
					BOOST_FOREACH( unsigned pre_s, predecessor_states[ s ] ) //for each predecessor state (or equivalently for every transition)
					{
						p_s_given_pre[ pre_s ] = p_r_given_predecessor[ pre_s ][ s ];
					}
					normalise_discrete_probabilities( p_s_given_pre );
					*/

					for( unsigned pre_s = 0; S != pre_s; ++pre_s )
					//BOOST_FOREACH( unsigned pre_s, predecessor_states[ s ] ) //for each predecessor state (or equivalently for every transition)
					{
						const transition_parameters & params = trans_params[ pre_s ][ s ];
						//VERIFY( params.get<0>() );

						const double p_pre_s = var_dist->rho[ n ][ i-1 ][ pre_s ];
						const bool g_pre_s = map.g( pre_s );
						const bool c_pre_s = map.c( pre_s );
						const unsigned k_pre_s = map.k( pre_s );
						const unsigned m_pre_s = map.m( pre_s );

						const double log_p_r_given_pre = params.get<0>() ? log_p_r_given_predecessor[ pre_s ][ s ] : -10.0 ;
						if( update_rho )
						{
							log_rho[ i ][ s ] += log_p_r_given_pre * p_pre_s;
							log_rho[ i-1 ][ pre_s ] += log_p_r_given_pre * p_s;
						}

						LL += log_p_r_given_pre * p_pre_s * p_s;
						VERIFY( ! MY_ISNAN( LL ) );

						if( params.get<0>() && update_tau )
						{
							//we don't want to update tau for sure thing transitions from gap bases
							if( ! map.g( pre_s ) )
							{
								//or from the last base
								if( ! (c_pre_s && 1 == k_pre_s) && ! (! c_pre_s && K == k_pre_s) )
								{
									const double a = params.get<1>();
									const double b = params.get<2>();
									const unsigned t_k = params.get<3>();
									var_dist->tau[ t_k ][ 0 ] += p_pre_s * p_s * a;
									var_dist->tau[ t_k ][ 1 ] += p_pre_s * p_s * ( a + b );
								}
							}
						}
					}
				}
			}
		}

		//update rho
		if( update_rho )
		{
			for( unsigned i = 0; I != i; ++i ) //for each base
			{
				probabilities_from_logs_of_choose( log_rho[ i ], var_dist->rho[ n ][ i ] );
			}
		}

	}

	//normalise eta if we are updating it
	if( update_eta )
	{
		for( unsigned m = 0; E != m; ++m ) //for each emission distribution
		{
			normalise_discrete_probabilities( var_dist->eta[ m ] );
		}
	}


	VERIFY( ! MY_ISNAN( LL ) );
	return LL;
}

void state_map_cached::init()
{
	const state_map_uncached map( K );
	const unsigned S = map.num_states();

	_g.resize( S );
	_c.resize( S );
	_b.resize( S );
	_k.resize( S );
	_m.resize( S );
	_s.resize( boost::extents[ 2 ][ 2 ][ K + 1 ] );

	for( unsigned s = 0; S != s; ++s )
	{
		const bool g = _g[ s ] = map.g( s );
		const bool c = _c[ s ] = map.c( s );
		_b[ s ] = map.b( s );
		const unsigned k = _k[ s ] = map.k( s );
		_m[ s ] = map.m( s );
		_s[ g ? 1 : 0 ][ c ? 1 : 0 ][ k ] = s;
	}
}

namespace impl
{
	typedef std::map< unsigned, state_map_cached > state_map_map;
	state_map_map state_maps;
}

const state_map & get_state_map( unsigned K )
{
	impl::state_map_map::iterator i = impl::state_maps.find( K );
	if( impl::state_maps.end() == i )
	{
		i = impl::state_maps.insert( std::make_pair( K, state_map_cached( K ) ) ).first;
	}
	return i->second;
}



} //namespace hmm
} //namespace gapped_pssm
} //namespace biopsy

