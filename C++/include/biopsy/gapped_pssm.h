/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_GAPPED_PSSM_H_
#define BIOPSY_GAPPED_PSSM_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include <biopsy/defs.h>

#include <bio/compressed_int_array.h>

#include <boost/multi_array.hpp>

#include <vector>


/** Compressed dna functions. */
namespace biopsy {

typedef compressed_int_array< unsigned, 2 > compressed_dna_seq;
typedef std::vector< compressed_dna_seq > dna_seq_list;

void string_to_compressed_dna( const std::string & str, compressed_dna_seq & seq );
void compressed_dna_to_string( const compressed_dna_seq & seq, std::string & str );

} //namespace biopsy



/** Dna as vector of unsigned chars. 0 is 'a', 1 is 'c', 2 is 'g', 3 is 't', 4 is 'n'. */
namespace biopsy {

typedef std::vector< unsigned char > dna_vec;
typedef std::vector< dna_vec > dna_vec_list;

void string_to_dna_vec( const std::string & str, dna_vec & seq );
void dna_vec_to_string( const dna_vec & seq, std::string & str );

} //namespace biopsy





namespace biopsy {
namespace gapped_pssm {



struct variational_model
{
	//typedefs
	typedef biopsy::dna_vec sequence;
	typedef biopsy::dna_vec_list seq_vec;
	typedef std::vector< double > vector;
	typedef std::vector< vector > vector_vec;
	typedef boost::multi_array< double, 2 > array;

	//members
	seq_vec seqs; /**< The sequences. */
	unsigned K; /**< Length of pssm (without gap). */
	vector alpha; /**< Beta prior parameters for gamma: the likelihood of a gap. */
	vector varphi; /**< Dirichlet prior parameters for background distribution. */
	vector phi; /**< Dirichlet prior parameters for pssm distribution. */
	vector lambda; /**< Variational parameter for gamma. */
	vector eta; /**< Variational parameter for location of the gap. */
	vector mu; /**< Variational parameter for g: has_gap variable. */
	array omega; /**< Variational parameters for background and pss distributions. */
	array log_p_x_given_r;
	array p_x_given_r;
	vector_vec nu; /**< Variational parameters for start positions of sites. */

	//methods
	variational_model( 
		unsigned K,
		const seq_vec & seqs,
		const vector & alpha,
		const vector & varphi,
		const vector & phi );

	void update();
	double log_likelihood() const;
	void shift( int offset );
	unsigned blank_sites( double l );
	void initialise_variational_params();

protected:
	double expectation_log_gamma() const;
	double expectation_log_1_minus_gamma() const;
	void recalc_after_omega_changes();
	void update_mu_nu_eta();
	void update_lambda();
	void update_omega();

	template< typename fn >
	fn
	generate_combinations( 
		unsigned n,
		fn f,
		bool only_sites = true ) const;
};


template< typename fn >
fn
variational_model::generate_combinations( 
	unsigned n,
	fn f,
	bool only_sites ) const
{
	const double p_g[2] = { 1.0 - mu[n], mu[n] };
	for( unsigned s = 0; nu[n].size() != s; ++s )
	{
		const unsigned s_position = s >> 1;
		const bool is_rev_comp = s & 1;
		const unsigned i_start = only_sites ? s_position : 0;
		const unsigned i_end = only_sites ? s_position + K + 1 : seqs[n].size();
		for( unsigned i = i_start; i_end != i; ++i )
		{
			const unsigned base_at_i = seqs[n][i];
			const unsigned x = ( 4 == base_at_i ? 4 : (is_rev_comp ? 3 - base_at_i : base_at_i) );
			const bool outside_pssm = i < s_position || i > s_position + K;

			for( unsigned h = 0; eta.size() != h; ++h )
			{
				for( unsigned g = 0; 2 != g; ++g )
				{
					// calculate r - index into pssm - 0 is background
					unsigned r;
					if( outside_pssm )
					{
						r = 0;
					} 
					else
					{
						const unsigned site_idx = 
							is_rev_comp
								? K + s_position - i
								: i - s_position;
						const unsigned gap_position = g ? h : K - 1;
						if( site_idx == gap_position + 1 ) r = 0;
						else if( site_idx > gap_position ) r = site_idx;
						else r = site_idx + 1;
						assert( 0 <= r && r < K + 1 );
					}

					//ignore anywhere outside the pssm bases
					if( only_sites && 0 == r )
					{
						continue;
					}

					f( 
						s,
						nu[n][s],
						h,
						eta[h],
						0 == g ? false : true,
						p_g[g],
						i,
						x,
						r );
				}
			}
		}
	}
	return f;
}


} //namespace gapped_pssm
} //namespace biopsy

#endif //BIOPSY_GAPPED_PSSM_H_
