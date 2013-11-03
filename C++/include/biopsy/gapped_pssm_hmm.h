/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_GAPPED_PSSM_HMM_H_
#define BIOPSY_GAPPED_PSSM_HMM_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/gsl.h"
#include "biopsy/containers.h"

#include <boost/dynamic_bitset.hpp>

namespace biopsy {
namespace gapped_pssm {

typedef unsigned char dna_base;
typedef std::vector< dna_base > sequence;
typedef std::vector< sequence > sequence_vec;



struct observed_sequences
{
	typedef boost::shared_ptr< observed_sequences > ptr;

	sequence_vec _sequences;
	unsigned N() const { return _sequences.size(); }
	unsigned I( unsigned n ) const { return _sequences[n].size(); }
	unsigned x( unsigned n, unsigned i ) const { return _sequences[n][i]; }
};



namespace hmm {


/** Maps states to g,c,b,k and back again. */
struct state_map_uncached
{
	unsigned K;

	state_map_uncached( unsigned K ) : K( K ) { }

	inline unsigned num_states() const { return 4 * K - 1; }

	inline bool g( unsigned s ) const { 
		if( 0 == s ) return false;
		return bool(((s - 1) % (2 * K - 1)) % 2);
	}

	inline bool c( unsigned s ) const { return s >= 2 * K; } //is s a reverse complement state?

	inline bool b( unsigned s ) const { return s == 0; } //is s the background state?

	inline unsigned k( unsigned s ) const { 
		if( b(s) ) return 0;
		else return ( (s - 1) % (2 * K - 1) ) / 2 + 1;
	}

	inline unsigned m( unsigned s ) const {
		if( b(s) ) return 0;
		return 2 * k(s) - ( g(s) ? 0 : 1 );
	}

	inline unsigned s( bool g, bool c, unsigned k ) const {
		if( 0 == k ) return 0;
		k = 2*(k-1)+1;
		if( g ) ++k;
		if( c ) k += 2 * K - 1;
		return k;
	}
};


/** Maps states to g,c,b,k and back again. */
struct state_map_cached
{
	typedef boost::dynamic_bitset<> dynamic_bitset;

	unsigned K;
	dynamic_bitset _g;
	dynamic_bitset _c;
	dynamic_bitset _b;
	unsigned_vector _k;
	unsigned_vector _m;
	unsigned_array_3d _s;

	inline unsigned num_states() const { return 4 * K - 1; }
	inline bool g( unsigned s ) const { return _g[ s ]; }
	inline bool c( unsigned s ) const { return _c[ s ]; } //is s a reverse complement state?
	inline bool b( unsigned s ) const { return _b[ s ]; } //is s the background state?
	inline unsigned k( unsigned s ) const { return _k[ s ]; }
	inline unsigned m( unsigned s ) const { return _m[ s ]; }
	inline unsigned s( bool g, bool c, unsigned k ) const { return _s[ g ? 1 : 0 ][ c ? 1 : 0 ][ k ]; }

	state_map_cached( unsigned K ) : K( K ) 
	{
		init();
	}

protected:
	void init();
};

typedef state_map_cached state_map;
const state_map & get_state_map( unsigned K );

void calculate_predecessor_states( unsigned K, unsigned_vector_vec & predecessor_states );


struct observed_data
{
	typedef boost::shared_ptr< observed_data > ptr;

	observed_sequences X;
	unsigned K;
	double_vector Psi;
	double_vector Theta;
	boost::array< double, 2 > Phi;
	boost::array< double, 2 > Upsilon;
	
	observed_data(
		unsigned K,
		const observed_sequences & X,
		double_vector Psi,
		double_vector Theta,
		double_vector Phi,
		double_vector Upsilon );

	inline unsigned S() const { return get_state_map( K ).num_states(); }
	inline unsigned E() const { return 2 * K; }
};


/** A transition probability is a function of t := a0 + a1*t[k] */
typedef boost::tuple< bool /* has transition */, double /* a0 */, double /* a1 */, unsigned /* k */ > transition_parameters;
typedef boost::multi_array< transition_parameters, 2 > transition_parameter_array;

struct hidden_data
{
	double_array e; /**< emission probs */
	double_vector t; /**< transition probs */

	/** Make a draw from the distribution defined by the observed_data (ignoring the sequences). */
	hidden_data( const observed_data & data );

	/** Draw sequence of r's of length I from distribution defined by e and t. */
	void draw_r( unsigned I, unsigned_vector & r ) const;

	/** Draw sequence of x's of length I from distribution defined by e and t. */
	void draw_sequence( unsigned I, unsigned_vector & r, sequence & seq ) const;

protected:
	transition_parameter_array trans_params;
};


struct variational_distribution
{
	typedef boost::shared_ptr< variational_distribution > ptr;

	double_vector_vec_vec rho; /**< variational dist over hidden variable r : the state for each base. */
	double_array eta; /**< variational dist over hidden variable e : the emission probabilities. */
	double_array tau; /**< variational dist over hidden variable t : the transition probabilities. */

	variational_distribution( const observed_data & data );

	inline double expected_t( unsigned k ) const { return tau[k][1] / ( tau[k][0] + tau[k][1] ); }
	void r_mode( unsigned_vector_vec & r ) const;
};

struct model
{
	observed_data::ptr data;
	variational_distribution::ptr var_dist;
	unsigned_vector_vec predecessor_states; /**< Those states that can preceed any given state. */
	transition_parameter_array trans_params;


	model( observed_data::ptr data );

	/** Update the variational distribution parameters. Returns log likelihood. */
	double update(
		bool update_rho = true, 
		bool update_eta = true, 
		bool update_tau = true );
};


void calculate_p_r_given_predecessor(
	const transition_parameter_array & trans_params,
	const variational_distribution & var_dist,
	double_array & p_r_given_predecessor );

void calculate_logs(
	const double_array & p,
	double_array & log_p );

void calculate_p_predecessor_given_r( 
	const double_array & p_r_given_predecessor,
	double_array & p_predecessor_given_r );

} //namespace hmm
} //namespace gapped_pssm_hmm
} //namespace biopsy

#endif //BIOPSY_GAPPED_PSSM_HMM_H_
