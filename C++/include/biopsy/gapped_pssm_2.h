/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_GAPPED_PSSM_2_H_
#define BIOPSY_GAPPED_PSSM_2_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/gsl.h"
#include "biopsy/gapped_pssm.h"
#include "biopsy/containers.h"







namespace biopsy {
namespace gapped_pssm_2 {

typedef biopsy::dna_vec sequence; /**< A dna sequence. */
typedef biopsy::dna_vec_list seq_vec; /**< A list of dna sequences. */


struct observed_data
{
	typedef boost::shared_ptr< observed_data > ptr;

	dna_vec_list		X;				/**< The observed sequences. */
	double_vector		A;				/**< Prior for s : does site start at given base? */
	double_vector		B;				/**< Prior for t : does given site have gap? */
	double_vector		V;				/**< Prior for part of w : the gap base distribution. */
	double_vector		W;				/**< Prior for part of w : the pssm base distribution. */
	unsigned			K;				/**< Length of pssm (excluding gap). */

	/** Construct from observed sequences and hyper-parameters. */
	observed_data( 
		const dna_vec_list & sequences,
		unsigned K,
		double expected_num_sites_per_seq,
		double strength_of_belief_in_num_sites_per_seq,
		double strength_of_gap_dist,
		double strength_of_pssm_dist );

	unsigned N() const; /** The number of sequences. */
	unsigned L( unsigned n ) const; /** The length of sequence n. */

};

struct variational_distribution
{
	//members
	double_vector_vec	alpha;			/**< Variational distribution over a : does site start at base i in sequence n? */
	double_vector_vec	beta;			/**< Variational distribution over b : does site at base i in sequence n have gap? */
	double_vector_vec	gamma;			/**< Variational distribution over c : is site at base i in sequence n a reverse complement? */
	double_vector		eta;			/**< Variational distribution over e : where the gap is. */
	double_vector		tau;			/**< Variational distribution over t : prob. any given site has gap. */
	double_array		omega;			/**< Variational distribution over w : the base distribution. */

	variational_distribution(
		const observed_data & data );

	double expected_t() const;							/**< Expected value of t under variational distribution. */
	void calc_p_e( double_vector & p_w ) const;			/**< p(e) under variational distribution. */
	void calc_expected_w( double_array & w ) const;		/**< p(e) under variational distribution. */
};

struct variational_model
{
	observed_data				data;					/**< The observed data. */
	variational_distribution	var_dist;				/**< The variational distribution. */
	double_array				p_x_given_r;			/**< p(x|r) cached for efficiency. */
	double_array				log_p_x_given_r;		/**< log p(x|r) cached for efficiency. */

	variational_model( const observed_data & data );

	/** Iterate through each combination of a, b, c and g. */
	template< typename fn >
	fn
	generate_combinations( 
		unsigned n,
		fn f ) const;

	double update();							/**< Update the variational distribution parameters. Returns log likelihood. */

protected:
	void update_p_x_given_r_expectations();
};

inline
unsigned calc_r(
	unsigned K,
	unsigned gap_position,
	unsigned offset,
	bool has_gap,
	bool rev_comp )
{
	//reverse everything if rev_comp
	if( rev_comp ) offset = K - offset;

	//if no gap - set gap position to end
	if( ! has_gap ) gap_position = K - 1;

	//where are relative to adjusted gap position?
	if( offset <= gap_position ) return offset + 1; //below gap
	if( offset == gap_position + 1 ) return 0; //0 means are in gap
	return offset; //above gap
};

inline
unsigned calc_offset(
	unsigned K,
	unsigned gap_position,
	unsigned r,
	bool has_gap,
	bool rev_comp )
{
	//if no gap - set gap position to end
	if( ! has_gap ) gap_position = K - 1;

	//where are we relative to the gap?
	if( 0 == r ) r = gap_position + 1; //we're in the gap
	else if( r < gap_position + 1 ) --r;

	//reverse everything if rev_comp
	if( rev_comp ) return K - r;
	else return r;
};

template< typename fn >
fn
variational_model::generate_combinations( 
	unsigned n,
	fn f ) const
{
	double_vector p_e_vec;
	var_dist.calc_p_e( p_e_vec );

	//for each putative site start position
	if( data.L( n ) < data.K + 1 ) throw std::logic_error( "Sequence too short for pssm" );
	const unsigned start_position_end = data.L( n ) - data.K - 1;
	for( unsigned a = 0; start_position_end != a; ++a )
	{
		const double p_a = var_dist.alpha[n][a];
		const unsigned site_end = a + data.K + 1;

		//for each instantiation with and without gap
		for( unsigned b = 0; 2 != b; ++b )
		{
			const double p_b = b ? var_dist.beta[n][a] : ( 1.0 - var_dist.beta[n][a] );

			//for each gap position
			for( unsigned e = 0; data.K != e; ++e )
			{
				const double p_e = p_e_vec[ e ];

				//for each orientation - reverse complement or not
				for( unsigned c = 0; 2 != c; ++c )
				{
					const double p_c = c ? var_dist.gamma[n][a] : ( 1.0 - var_dist.gamma[n][a] );

					//for each base in putative site
					for( unsigned i = a; site_end != i; ++i )
					{
						//index into pssm
						const unsigned r = calc_r( data.K, e, i - a, bool( b ), bool( c ) );

						//is the base complemented?
						const unsigned x =
							c
								? 3 - data.X[ n ][ i ] //yes
								: data.X[ n ][ i ] //no
							;

						f( i, a, p_a, b, p_b, c, p_c, e, p_e, r, x );
					}
				}
			}
		}
	}
	return f;
}

} //namespace gapped_pssm
} //namespace biopsy

#endif //BIOPSY_GAPPED_PSSM_2_H_
