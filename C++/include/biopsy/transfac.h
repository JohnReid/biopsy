/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_TRANSFAC_H_
#define BIOPSY_TRANSFAC_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"
#include "biopsy/pssm.h"



namespace bio {
struct BiobasePssmFilter;
struct BiobaseTablePssmEntry;
}


namespace biopsy
{

/** Is it a TRANSFAC PSSM name? */
bool is_transfac_pssm( const std::string & pssm_name );

/**
Get the biobase PSSM entry for the given PSSM name.
*/
bio::BiobaseTablePssmEntry * 
get_transfac_pssm_entry( const std::string & pssm_name );

/**
Get the biobase accession id for a given PSSM name
*/
std::string
get_transfac_pssm_accession( 
	const std::string & pssm_name );
	

/**
Get the name for a given PSSM
*/
std::string
get_transfac_pssm_name(
	const std::string & pssm );
	

/**
Get the sequences a PSSM was built from.
*/
string_vec_ptr
get_transfac_pssm_sequences(
	const std::string & pssm );


/**
Get the minimise false positives threshold for a PSSM.
*/
double
get_transfac_pssm_min_fp_threshold(
	const std::string & pssm );


/**
Get the minimise false negatives threshold for a PSSM.
*/
double
get_transfac_pssm_min_fn_threshold(
	const std::string & pssm );


/**
Get the minimise sum false positives and negatives threshold for a PSSM.
*/
double
get_transfac_pssm_min_sum_threshold(
	const std::string & pssm );


/**
Get the TRANSFAC PSSMs' accessions.
*/
string_vec_ptr
get_transfac_pssm_accessions( const bio::BiobasePssmFilter & filter );


/**
Get the default filter for TRANSFAC PSSMs. I.e. use consensus sequences and only vertebrate matrices.
*/
bio::BiobasePssmFilter
get_default_transfac_pssm_filter();


/**
Get the factors linked to the PSSM.
*/
string_vec_ptr
get_factors_for_pssm( const std::string & pssm_acc );


/**
Get the PSSMs for the factor.
*/
string_vec_ptr
get_pssms_for_factor( const std::string & factor_acc );


/**
Get the matrices in TRANSFAC
*/
string_vec_ptr get_transfac_matrices( );

/**
Get the sites in TRANSFAC
*/
string_vec_ptr get_transfac_sites( );

/**
Get the factors in TRANSFAC
*/
string_vec_ptr get_transfac_factors( );

/**
Get the fragments in TRANSFAC
*/
string_vec_ptr get_transfac_fragments( );

/**
Get the genes in TRANSFAC
*/
string_vec_ptr get_transfac_genes( );

/**
Get the factors for the fragment.
*/
string_vec_ptr
get_factors_for_fragment( const std::string & fragment_acc );

/**
Get the fragments for the factor
*/
string_vec_ptr
get_fragments_for_factor( const std::string & factor_acc );

/**
Get the fragment's sequence
*/
std::string
get_fragment_sequence( const std::string & fragment_acc );


} //namespace biopsy

#endif //BIOPSY_TRANSFAC_H_
