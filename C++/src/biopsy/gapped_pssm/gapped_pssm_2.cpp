/**
@file

Copyright John Reid 2006

*/

#define GSL_DLL

#include <biopsy/gapped_pssm_2.h>
#include <biopsy/gsl.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>



namespace biopsy {
namespace gapped_pssm_2 {

observed_data::observed_data( 
	const dna_vec_list & sequences,
	unsigned K,
	double expected_num_sites_per_seq,
	double strength_of_belief_in_num_sites_per_seq,
	double strength_of_gap_dist,
	double strength_of_pssm_dist )
: X( sequences )
, A( sequences.size() )
, B( 2, 1.0 )
, V( 4, strength_of_gap_dist )
, W( 4, strength_of_pssm_dist )
, K( K )
{
	for( unsigned n = 0; sequences.size() != n; ++n )
	{
		if( K + 1 >= X[n].size() ) throw std::logic_error( BIOPSY_MAKE_STRING( "Gapped pssm longer than sequence " << n ) );
		A[ n ] = expected_num_sites_per_seq * strength_of_belief_in_num_sites_per_seq / ( X[n].size() - K - 1 );
	}
}

unsigned
observed_data::N() const
{
	return X.size();
}

unsigned
observed_data::L( unsigned n ) const
{
	return X[n].size();
}

namespace impl {

template< 
	typename log_vector, 
	typename out_vector 
>
void
probabilities_from_logs_of_choose( 
	const log_vector & logs,
	out_vector & out )
{
	// make the largest 0 - as we can get underflow problems
	const double largest = *( std::max_element( logs.begin(), logs.end() ) );
	out.resize( logs.size() );
	for( unsigned i = 0; logs.size() != i; ++i )
	{
		const double v = logs[i] - largest;
		out[i] = (v <= GSL_LOG_DBL_MIN) ? 0.0 : gsl_sf_exp( v );
	}
	const double total = std::accumulate( out.begin(), out.end(), 0.0 );
	for( unsigned i = 0; out.size() != i; ++i )
	{
		out[i] /= total;
	}
}

void
initialise_bernoulli_log_array_empty( 
	double_array & a )
{
	for( unsigned i = 0; a.size() != i; ++i )
	{
		a[i][1] = 0.0;
		a[i][0] = 0.0;
	}
}

void
initialise_bernoulli_log_array( 
	double_array & a,
	double p )
{
	for( unsigned i = 0; a.size() != i; ++i )
	{
		a[i][1] = gsl_sf_log( p );
		a[i][0] = gsl_sf_log( 1.0 - p );
	}
}

template< 
	typename log_vector, 
	typename out_vector 
>
void
probabilities_from_logs_of_bernoulli( 
	const log_vector & logs,
	out_vector & out )
{
	out.resize( logs.size() );
	for( unsigned i = 0; logs.size() != i; ++i )
	{
		double l_true = logs[ i ][ 1 ];
		double l_false = logs[ i ][ 1 ];
		const double largest = std::max( l_true, l_false ); //avoid underflow
		l_true -= largest;
		l_false -= largest;
		const double p_true = ( l_true <= GSL_LOG_DBL_MIN ) ? 0.0 : gsl_sf_exp( l_true );
		const double p_false = ( l_false <= GSL_LOG_DBL_MIN ) ? 0.0 : gsl_sf_exp( l_false );
		out[i] = p_true / ( p_false + p_true );
	}
}

} //namespace impl

variational_distribution::variational_distribution(
	const observed_data & data )
: alpha( data.N() )
, beta( data.N() )
, gamma( data.N() )
, eta( data.K, 1.0 / data.K )
, tau( 2, 0.5 )
, omega( boost::extents[ data.K + 1 ][ 4 ] )
{
	for( unsigned n = 0; data.N() != n; ++n )
	{
		if( data.L(n) <= data.K + 1 ) throw std::logic_error( "Sequence too small for pssm" );
		
		const unsigned num_possible_site_start_positions = data.L(n) - data.K - 1;
		
		alpha[n].resize( num_possible_site_start_positions );
		beta[n].resize( num_possible_site_start_positions );
		gamma[n].resize( num_possible_site_start_positions );

		for( unsigned i = 0; num_possible_site_start_positions != i; ++i )
		{
			alpha[n][i] = gsl_ran_beta(
				get_gsl_rng(),
				1,
				num_possible_site_start_positions - 1 );
			beta[n][i] = gsl_ran_beta(
				get_gsl_rng(),
				1,
				1 );
			gamma[n][i] = gsl_ran_beta(
				get_gsl_rng(),
				1,
				1 );
		}

	}

	//initialise omega
	omega[0][0] = omega[0][1] = omega[0][2] = omega[0][3] = 100.0;
	for( unsigned j = 1; data.K + 1 != j; ++j )
	{
		omega[j][0] = omega[j][1] = omega[j][2] = omega[j][3] = 0.01;
	}
}

double
variational_distribution::expected_t() const
{
	return tau[1] / ( tau[0] + tau[1] );
}

void
variational_distribution::calc_expected_w( double_array & w ) const
{
	w.resize( boost::extents[ omega.size() ][ 4 ] );
	for( unsigned j = 0; omega.size() != j; ++j )
	{
		const double sum = std::accumulate( omega[ j ].begin(), omega[ j ].end(), 0.0 );
		for( unsigned x = 0; 4 != x; ++x )
		{
			w[ j ][ x ] = omega[ j ][ x ] / sum;
		}
	}
}

void
variational_distribution::calc_p_e( double_vector & p_e ) const
{
	p_e.resize( eta.size() );
	const double sum = std::accumulate( eta.begin(), eta.end(), 0.0 );
	for( unsigned e = 0; p_e.size() != e; ++e )
	{
		p_e[ e ] = eta[ e ] / sum;
	}
}


variational_model::variational_model( const observed_data & data )
: data( data )
, var_dist( data )
, p_x_given_r( boost::extents[ data.K + 1 ][ 4 ] )
, log_p_x_given_r( boost::extents[ data.K + 1 ][ 5 ] )
{
	update_p_x_given_r_expectations();
}

void
variational_model::update_p_x_given_r_expectations()
{
	for( unsigned r = 0; data.K + 1 != r; ++r )
	{
		double_array::reference lp = log_p_x_given_r[ r ];
		const double t = std::accumulate( var_dist.omega[r].begin(), var_dist.omega[r].end(), 0.0 );
		const double psi_t = gsl_sf_psi( t );
		for( unsigned x = 0; 4 != x; ++x )
		{
			lp[x] = gsl_sf_psi( var_dist.omega[r][x] ) - psi_t;
			p_x_given_r[r][x] = var_dist.omega[r][x] / t;
		}
		lp[4] = std::accumulate( 
			lp.begin(),
			lp.begin() + 4,
			0.0 ) * 0.25;
	}
}

struct update_fn
{
	double						LL;
	double_array &				log_alpha;
	double_array &				log_beta;
	double_array &				log_gamma;
	double_vector &				log_eta;
	double_array &				omega;
	const double_array &		log_p_x_given_r;

	update_fn(
		const double_array & log_p_x_given_r,
		double_array & log_alpha,
		double_array & log_beta,
		double_array & log_gamma,
		double_vector & log_eta,
		double_array & omega )
		: LL( 0.0 )
		, log_alpha( log_alpha )
		, log_beta( log_beta )
		, log_gamma( log_gamma )
		, log_eta( log_eta )
		, omega( omega )
		, log_p_x_given_r( log_p_x_given_r )
	{
	}

	inline
	void 
	operator()(
		unsigned i, 
		unsigned a,
		double p_a,
		unsigned b,
		double p_b,
		unsigned c,
		double p_c,
		unsigned e,
		double p_e,
		unsigned r,
		unsigned x )
	{
		//weights in expectation
		const double alpha_weight = p_b * p_c * p_e;
		const double beta_weight = p_a * p_c * p_e;	
		const double gamma_weight = p_a * p_b * p_e;
		const double eta_weight = p_a * p_b * p_c;	
		const double omega_weight = p_a * p_b * p_c * p_e;

		const double expectation = log_p_x_given_r[ r ][ x ];

		//update alpha, beta, gamma, eta
		//std::cout << i << "," << a << "," << b << "," << c << "\n";
		log_alpha[ a ][ 1 ] += alpha_weight * expectation;
		log_beta[ a ][ b ] += beta_weight * expectation;
		log_gamma[ a ][ c ] += gamma_weight * expectation;
		log_eta[ e ] += eta_weight * expectation;

		//update omega
		if( 4 == x ) //x is an 'N'
		{
			for( unsigned j = 0; 4 != j; ++j )
			{
				omega[r][j] += omega_weight * 0.25;
			}
		}
		else
		{
			omega[r][x] += omega_weight;
		}
	}
};



double
variational_model::update()
{
	static const double log_quarter = gsl_sf_log( 0.25 );

	double LL = 0.0;

	//initialise log eta array
	double_vector log_eta( var_dist.eta.size(), 0.0 );

	//for each sequence
	for( unsigned n = 0; data.N() != n; ++n )
	{
		//initialise log alpha array
		double_array log_alpha( boost::extents[ var_dist.alpha[ n ].size() ][ 2 ] );
		impl::initialise_bernoulli_log_array( log_alpha, data.A[n] );
		for( unsigned a = 0; log_alpha.size() != a; ++a ) log_alpha[ a ][ 0 ] = log_quarter;

		//initialise log beta array
		double_array log_beta( boost::extents[ var_dist.beta[ n ].size() ][ 2 ] );
		impl::initialise_bernoulli_log_array( log_beta, var_dist.expected_t() );

		//initialise log gamma array
		double_array log_gamma( boost::extents[ var_dist.gamma[ n ].size() ][ 2 ] );
		impl::initialise_bernoulli_log_array_empty( log_gamma );

		//initialise omega to update
		for( unsigned j = 0; data.K + 1 != j; ++j )
		{
			if( 0 == j ) for( unsigned x = 0; 4 != x; ++x ) var_dist.omega[j][x] = data.V[x]; //init from V
			else for( unsigned x = 0; 4 != x; ++x ) var_dist.omega[j][x] = data.W[x]; //init from W
		}

		//go through each combination
		LL += generate_combinations( 
			n,
			update_fn( 
				log_p_x_given_r, 
				log_alpha, 
				log_beta, 
				log_gamma, 
				log_eta, 
				var_dist.omega ) ).LL;

		//convert logs back to probabilities
		impl::probabilities_from_logs_of_bernoulli( log_alpha, var_dist.alpha[ n ] );
		impl::probabilities_from_logs_of_bernoulli( log_beta, var_dist.beta[ n ] );
		impl::probabilities_from_logs_of_bernoulli( log_gamma, var_dist.gamma[ n ] );
	}

	impl::probabilities_from_logs_of_choose( log_eta, var_dist.eta );

	//omega has changed so...
	update_p_x_given_r_expectations();

	return LL;
}

} //namespace gapped_pssm_2
} //namespace biopsy

