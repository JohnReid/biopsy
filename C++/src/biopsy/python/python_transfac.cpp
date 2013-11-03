/**
@file

Copyright John Reid 2006

*/

#include <boost/python.hpp>
#include "biopsy/transfac.h"
#include "biopsy/python.h"

#include <bio/biobase_filter.h>


using namespace boost::python;
using namespace std;


namespace biopsy {




//
// export
//
void export_transfac()
{
#if 0
	class_< BIO_NS::BiobasePssmFilter >( "TransfacPssmFilter", boost::python::init< bool, const std::string &, const std::string & >() )
		.def(
			"all_pssms", 
			&BIO_NS::BiobasePssmFilter::get_all_pssms_filter )
		.staticmethod( "all_pssms" )
		;
#endif

	def(
		"get_default_transfac_pssm_filter", 
		get_default_transfac_pssm_filter,
		"The default filter to limit which pssms from transfac we use" );
	def(
		"get_transfac_pssm_accession", 
		get_transfac_pssm_accession,
		"Get the pssms accession from its name" );
	def( 
		"get_transfac_pssm_name",
		get_transfac_pssm_name,
		"Get the transfac pssm's name from its accession" );
	def( "get_transfac_pssm_sequences", get_transfac_pssm_sequences );
	def(
		"get_transfac_pssm_accessions", 
		get_transfac_pssm_accessions,
		"Get the transfac accessions for the given filter" );
	def( "get_transfac_pssm_min_fp_threshold", get_transfac_pssm_min_fp_threshold );
	def( "get_transfac_pssm_min_fn_threshold", get_transfac_pssm_min_fn_threshold );
	def( "get_transfac_pssm_min_sum_threshold", get_transfac_pssm_min_sum_threshold );
	def( "get_factors_for_pssm", get_factors_for_pssm );
	def( "get_pssms_for_factor", get_pssms_for_factor );
	def( 
		"get_transfac_matrices", 
		get_transfac_matrices,
		"All the matrices in transfac" );
	def( 
		"get_transfac_sites", 
		get_transfac_sites,
		"All the sites in transfac" );
	def( 
		"get_transfac_factors", 
		get_transfac_factors,
		"All the factors in transfac" );
	def( 
		"get_transfac_fragments", 
		get_transfac_fragments,
		"All the fragments in transfac" );
	def( 
		"get_transfac_genes", 
		get_transfac_genes,
		"All the genes in transfac" );
	def( 
		"get_factors_for_fragment", 
		get_factors_for_fragment,
		"Get the factors for the given factor" );
	def( 
		"get_fragments_for_factor", 
		get_fragments_for_factor,
		"Get the fragments which contain the factor" );
	def( 
		"get_fragment_sequence", 
		get_fragment_sequence,
		"Get the sequence for the factor" );



}




} //namespace biopsy


