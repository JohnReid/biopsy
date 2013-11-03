/**
@file

Copyright John Reid 2006

*/

#include <biopsy/gapped_pssm.h>
#include <biopsy/gsl.h>

#include <vector>
#include <sstream>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>


namespace biopsy {


void string_to_compressed_dna( const std::string & str, compressed_dna_seq & seq )
{
	seq.resize( str.size() );
	for( unsigned i = 0; str.size() != i; ++i )
	{
		switch( str[i] )
		{
		case 'a':
		case 'A':
			seq.set( i, 0 );
			break;
		case 'c':
		case 'C':
			seq.set( i, 1 );
			break;
		case 'g':
		case 'G':
			seq.set( i, 2 );
			break;
		case 't':
		case 'T':
			seq.set( i, 3 );
			break;
		default:
			throw std::logic_error( "Non-nucleotide character in sequence" );
		}
	}
}

void compressed_dna_to_string( const compressed_dna_seq & seq, std::string & str )
{
	str.resize( seq.size );
	for( unsigned i = 0; str.size() != i; ++i )
	{
		switch( seq.get< unsigned >( i ) )
		{
		case 0: str[i] = 'a'; break;
		case 1: str[i] = 'c'; break;
		case 2: str[i] = 'g'; break;
		case 3: str[i] = 't'; break;
		default:
			throw std::logic_error( "Impossible value in compressed dna" );
		}
	}
}

} //namespace biopsy





namespace biopsy {

void dna_vec_to_string( const dna_vec & seq, std::string & str )
{
	str.resize( seq.size() );
	for( unsigned i = 0; str.size() != i; ++i )
	{
		switch( seq[i] )
		{
		case 0: str[i] = 'a'; break;
		case 1: str[i] = 'c'; break;
		case 2: str[i] = 'g'; break;
		case 3: str[i] = 't'; break;
		case 4: str[i] = 'n'; break;
		default:
			throw std::logic_error( "Impossible value in dna vec" );
		}
	}
}


void string_to_dna_vec( const std::string & str, dna_vec & seq )
{
	seq.resize( str.size() );
	for( unsigned i = 0; str.size() != i; ++i )
	{
		switch( str[i] )
		{
		case 'a': case 'A': seq[i] = 0; break;
		case 'c': case 'C': seq[i] = 1; break;
		case 'g': case 'G': seq[i] = 2; break;
		case 't': case 'T': seq[i] = 3; break;
		case 'n': case 'N': seq[i] = 4; break;
		default:
			throw std::logic_error( "Non-nucleotide character in sequence" );
		}
	}
}

} //namespace biopsy




namespace biopsy {
namespace gapped_pssm {


variational_model::variational_model(
	unsigned K,
	const seq_vec & seqs,
	const vector & alpha,
	const vector & varphi,
	const vector & phi )
	: seqs( seqs )
	, K( K )
	, alpha( alpha )
	, varphi( varphi )
	, phi( phi )
	, lambda( 2, 1.0 )
	, eta( K - 1, 1.0 / double( K - 1 ) )
	, mu( seqs.size(), 1.0 / double( seqs.size() ) )
	, omega( boost::extents[ K+1 ][ 4 ] )
	, log_p_x_given_r( boost::extents[ K+1 ][ 5 ] )
	, p_x_given_r( boost::extents[ K+1 ][ 4 ] )
	, nu( seqs.size() )
{
	initialise_variational_params();
}

void
variational_model::initialise_variational_params()
{
	//initialise variational parameters over start position
	for( unsigned i = 0; seqs.size() != i; ++i )
	{
		if( seqs[i].size() <= K ) throw std::logic_error( "Sequence too short for pssm" );
		const unsigned nu_size = 2 * ( seqs[i].size() - K ); //need 2* for reverse complement
		nu[i].resize( nu_size );
		std::fill( nu[i].begin(), nu[i].end(), 1.0 / double( nu_size ) );
	}

	//initialise var. params over background and pssm
	for( unsigned r = 0; K + 1 != r; ++r )
	{
		if( 0 == r ) 
		{
			std::fill( omega[r].begin(), omega[r].end(), 0.25 );
		}
		else
		{
			gsl_ran_dirichlet( 
				get_gsl_rng(), 
				omega[r].size(),
				&(*phi.begin()),
				&(*omega[r].begin()) );
		}
	}
	recalc_after_omega_changes();
}

void
variational_model::update()
{
	update_omega();
	update_mu_nu_eta();
	update_lambda();
}

double
variational_model::expectation_log_gamma() const
{
	return gsl_sf_psi( lambda[0] )
		- gsl_sf_psi( lambda[0] + lambda[1] );
}


double
variational_model::expectation_log_1_minus_gamma() const
{
	return gsl_sf_psi( lambda[1] )
		- gsl_sf_psi( lambda[0] + lambda[1] );
}


void
variational_model::recalc_after_omega_changes()
{
	for( unsigned r = 0; K + 1 != r; ++r )
	{
		array::reference lp = log_p_x_given_r[r];
		const double t = std::accumulate( omega[r].begin(), omega[r].end(), 0.0 );
		const double psi_t = gsl_sf_psi( t );
		for( unsigned x = 0; 4 != x; ++x )
		{
			lp[x] = gsl_sf_psi( omega[r][x] ) - psi_t;
			p_x_given_r[r][x] = omega[r][x] / t;
		}
		lp[4] = std::accumulate( 
			lp.begin(),
			lp.begin() + 4,
			0.0 ) * 0.25;
	}
}


template< typename log_vector, typename out_vector >
void
probabilities_from_logs( 
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


struct update_mu_nu_eta_fn
{
	typedef variational_model::vector vector;
	typedef variational_model::array array;

	vector & log_h;
	vector & log_s;
	vector & log_g;
	const array & log_p_x_given_r;

	update_mu_nu_eta_fn( 
		vector & log_h,
		vector & log_s,
		vector & log_g,
		const array & log_p_x_given_r )
		: log_h( log_h )
		, log_s( log_s )
		, log_g( log_g )
		, log_p_x_given_r( log_p_x_given_r )
	{ }

	inline
	void 
	operator()(
		unsigned s,
		double q_s,
		unsigned h,
		double q_h,
		bool g,
		double q_g,
		unsigned i,
		unsigned x,
		unsigned r )
	{
		const double g_weight = q_s * q_h; // weight in expectation for g
		const double s_weight = q_g * q_h; // weight in expectation for s
		const double h_weight = q_g * q_s; // weight in expectation for h
		const double expectation = log_p_x_given_r[r][x];
		log_g[g ? 1 : 0] += g_weight * expectation;
		log_s[s] += s_weight * expectation;
		log_h[h] += h_weight * expectation;
	}
};

void
variational_model::update_mu_nu_eta()
{
	vector log_h( eta.size(), 0.0 );
	for( unsigned n = 0; seqs.size() != n; ++n )
	{
		vector log_g( 2 );
		log_g[0] = expectation_log_gamma();
		log_g[1] = expectation_log_1_minus_gamma();
		vector log_s( nu[n].size(), 0.0 );
		generate_combinations(
			n,
			update_mu_nu_eta_fn(
				log_h,
				log_s,
				log_g,
				log_p_x_given_r ) );
			
		// update mu
		vector p_g( 2 );
		probabilities_from_logs( log_g, p_g );
		mu[n] = p_g[1];

		// update nu
		probabilities_from_logs( log_s, nu[n] );
	}

	// update eta
	probabilities_from_logs( log_h, eta );
}

void
variational_model::update_lambda()
{
	const double t = std::accumulate( mu.begin(), mu.end(), 0.0 );
	lambda[0] = alpha[0] + t;
	lambda[1] = alpha[1] + seqs.size() - t;
}


struct update_omega_fn
{
	typedef variational_model::array array;

	array & omega;

	update_omega_fn( 
		array & omega )
		: omega( omega )
	{ }

	inline
	void 
	operator()(
		unsigned s,
		double q_s,
		unsigned h,
		double q_h,
		bool g,
		double q_g,
		unsigned i,
		unsigned x,
		unsigned r )
	{
		const double q = q_s * q_h * q_g;
		if( 4 == x ) //x is an 'N'
		{
			for( unsigned i = 0; 4 != i; ++i )
			{
				omega[r][i] += q * 0.25;
			}
		}
		else
		{
			omega[r][x] += q;
		}
	}
};

void
variational_model::update_omega()
{
	for( unsigned j = 0; K + 1 != j; ++j )
	{
		for( unsigned x = 0; 4 != x; ++x )
		{
			omega[j][x] = (0 == j) ? varphi[x] : phi[x];
		}
	}
	for( unsigned n = 0; seqs.size() != n; ++n )
	{
		generate_combinations(
			n,
			update_omega_fn( omega ) );
	}
		
	// update cached values
	recalc_after_omega_changes();
}

struct log_likelihood_fn
{
	typedef variational_model::array array;
	typedef variational_model::vector vector;

	vector likelihoods;
	vector p_r_not_0;
	const array & p_x_given_r;

	log_likelihood_fn( unsigned n, const array & p_x_given_r )
		: likelihoods( n, 0.0 )
		, p_r_not_0( n, 0.0 )
		, p_x_given_r( p_x_given_r )
	{ }

	inline
	void 
	operator()(
		unsigned s,
		double q_s,
		unsigned h,
		double q_h,
		bool g,
		double q_g,
		unsigned i,
		unsigned x,
		unsigned r )
	{
		const double q = q_s * q_h * q_g;
		p_r_not_0[i] += q;
		if( 4 == x ) //x is an 'N'
		{
			likelihoods[i] += q * 0.25;
		}
		else
		{
			likelihoods[i] += q * p_x_given_r[r][x];
		}
	}

	double log_likelihood( const variational_model::sequence & seq ) const
	{
		double total = 0.0;
		for( unsigned i = 0; likelihoods.size() != i; ++i )
		{
			array::const_reference p_x_given_r_0 = p_x_given_r[0];
			total += 
				gsl_sf_log( 
					likelihoods[i] 
					+ ( 1.0 - p_r_not_0[i] ) * ( 4 == seq[i] ? 0.25 : p_x_given_r_0[seq[i]] ) );
		}
		return total;
	}
};

double 
variational_model::log_likelihood() const
{
	double total_ll = 0.0;
	for( unsigned n = 0; seqs.size() != n; ++n )
	{
		total_ll +=
			generate_combinations(
				n,
				log_likelihood_fn( seqs[n].size(), p_x_given_r ),
				true ).log_likelihood( seqs[n] );
	}
	return total_ll;
}


template< typename Array, typename InitElement >
void shift_array( 
	const Array & old_array, 
	Array & new_array,
	bool zero_special,
	int offset,
	const InitElement & init_element )
{
	const int lowest_idx = zero_special ? 1 : 0;

	//the rest of the indices do shift
	for( int new_idx = 0; new_array.size() != unsigned( new_idx ); ++new_idx )
	{
		if( zero_special && 0 == new_idx )
		{
			//zero-th index does not shift
			new_array[0] = old_array[0];
			continue;
		}

		const int old_idx = new_idx - offset;

		//do we index into old array?
		new_array[new_idx] = 
			( lowest_idx <= old_idx && old_idx < int( old_array.size() ) )
				? old_array[old_idx] //yes
				: init_element; //no
	}
}

template< typename vector >
boost::multi_array< double, 1 > multi_array_from_vec( const vector & v )
{
	boost::multi_array< double, 1 > result( boost::extents[ v.size() ] );
	for( unsigned i = 0; v.size() != i; ++i )
	{
		result[i] = v[i];
	}
	return result;
}

template< typename vector >
double
normalise( vector & v )
{
	const double s = std::accumulate( v.begin(), v.end(), 0.0 );
	for( unsigned i = 0; v.size() != i; ++i )
	{
		v[i] /= s;
	}
	return s;
}

template< typename vector >
double find_threshold( const vector & p, double l )
{
	vector sorted_p( p );
	std::sort( sorted_p.begin(), sorted_p.end() );
	double total = 0.0;
	unsigned i;
	for( i = 0; sorted_p.size() != i; ++i )
	{
		total += sorted_p[i];
		if( total > l ) break;
	}
	return sorted_p[i];
}

unsigned 
variational_model::blank_sites( double l )
{
	unsigned blanked = 0;

	//for each sequence
	for( unsigned n = 0; nu.size() != n; ++n )
	{
		const double threshold = find_threshold( nu[n], l );
		for( unsigned i = 0; nu[n].size() != i; ++i )
		{
			if( nu[n][i] >= threshold )
			{
				//blank site
				const unsigned s_position = i >> 1;
				for( unsigned j = s_position; s_position + K + 1 != j; ++j )
				{
					seqs[n][j] = 4;
				}

				++blanked;
			}
		}
	}

	return blanked;
}

void 
variational_model::shift( int offset )
{
	//omega - pssm
	array new_omega( omega );
	shift_array( omega, new_omega, true, offset, multi_array_from_vec( phi ) );
	omega = new_omega;
	recalc_after_omega_changes();

	//eta - where gap is
	vector new_eta( eta );
	shift_array( eta, new_eta, false, offset, 1.0 / double( eta.size() ) );
	eta = new_eta;
	normalise( eta );

	//nu - where starting position is - has to shift other direction from pssm
	for( unsigned i = 0; nu.size() != i; ++i )
	{
		vector new_nu( nu[i].size() );
		shift_array( nu[i], new_nu, false, - 2 * offset, 1.0 / double( nu[i].size() ) );
		nu[i] = new_nu;
		normalise( nu[i] );
	}

}


} //namespace gapped_pssm
} //namespace biopsy

