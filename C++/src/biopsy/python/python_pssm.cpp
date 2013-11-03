/**
@file

Copyright John Reid 2006

*/

#include <boost/python.hpp>
#include "biopsy/binding_hits.h"
#include "biopsy/analyse.h"
#include "biopsy/pssm.h"
#include "biopsy/custom_pssm.h"
#include "biopsy/python.h"

using namespace boost;
using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;

namespace boost { namespace python { namespace indexing {
template<>
struct value_traits<biopsy::nucleo_dist> : value_traits< int >
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
template<>
struct value_traits<biopsy::pssm_ptr> : value_traits< int >
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
} } }


namespace biopsy {


void export_pssm()
{
	/** A distribution over nucleotides. */
	class_< nucleo_dist >( "NucleoDist", boost::python::init< double, double, double, double >() )
		.def( "get", &nucleo_dist::get )
		.def( "get_total", &nucleo_dist::get_total )
		.def( "get_freq", &nucleo_dist::get_freq )
		;


	/**
	Likelihoods.
	*/
	if( ! converter::registry::query( type_id< likelihoods >() ) ) {
		class_< likelihoods, likelihoods_ptr, boost::noncopyable >( "Likelihoods" )
			.def( container_suite< likelihoods >() )
		;
	}
	register_ptr_to_python< likelihoods_ptr >();
	def( "get_likelihood_index", get_likelihood_index );
	def( "get_likelihood", get_likelihood );
	def( "get_odds_ratio", get_odds_ratio );
	def( "get_p_binding", get_p_binding );
	def( "get_p_binding_from_score", get_p_binding_from_score );
	def( "get_p_binding_using_p_value", get_p_binding_using_p_value );
	def( "get_odds_ratio_from_p_binding", get_odds_ratio_from_p_binding );


	/**
	A pssm.
	*/
	class_< pssm >( "Pssm" )
		.def( container_suite< pssm >() )
		;
	register_ptr_to_python< pssm_ptr >();

	/**
	A vector of pssms.
	*/
	class_< pssm_ptr_vec >( "PssmVec" )
		.def( container_suite< pssm_ptr_vec >() )
		;
	register_ptr_to_python< pssm_ptr_vec_ptr >();

	/**
	A pssm and its info.
	*/
	class_< pssm_info >(
		"PssmInfo",
		boost::python::init< const nucleo_dist::vec &, double, int, pssm_ptr, likelihoods_ptr, likelihoods_ptr >()
	)
		.def_readonly( "pssm", &pssm_info::_pssm )
		.def_readonly( "dists", &pssm_info::_dists )
		.def_readonly( "counts", &pssm_info::_counts )
		.def_readonly( "pseudo_count", &pssm_info::_pseudo_count )
		.def_readonly( "sites", &pssm_info::_number_of_sites )
		.def( "get_likelihoods", &pssm_info::get_dist )
		;
	def( "create_pssm", create_pssm );
	def( "calculate_likelihoods_under_pssm", calculate_likelihoods_under_pssm );
	def( "calculate_likelihoods_under_background", calculate_likelihoods_under_background );
	def( "add_pssm_to_cache", add_pssm_to_cache );
	def( "save_pssm_cache_state", save_pssm_cache_state );
	def( "clear_pssm_cache", clear_pssm_cache );
	def( "get_pssm", get_pssm, return_value_policy< return_by_value >() );
	def( "get_pssm_name", get_pssm_name );
	def( "get_pssm_url", get_pssm_url );
	def( "score_pssm", score_pssm );

	class_< pssm_parameters, noncopyable >( "PssmParameters", no_init )
		ADD_STATIC_PROPERTY(pssm_parameters,pseudo_counts)
		ADD_STATIC_PROPERTY(pssm_parameters,likelihoods_size)
		ADD_STATIC_PROPERTY(pssm_parameters,calculate_likelihoods_map_size)
		ADD_STATIC_PROPERTY(pssm_parameters,binding_background_odds_prior)
		ADD_STATIC_PROPERTY(pssm_parameters,use_cumulative_dists)
		ADD_STATIC_PROPERTY(pssm_parameters,use_p_value)
		ADD_STATIC_PROPERTY(pssm_parameters,use_score)
		ADD_STATIC_PROPERTY(pssm_parameters,avg_phylo_bayes)
		ADD_STATIC_PROPERTY(pssm_parameters,max_chain_num_boxes_limit)
		ADD_STATIC_PROPERTY(pssm_parameters,min_related_evidence_fraction)
		;


	/**
	A custom pssm and its info.
	*/
	class_< custom_pssm >( "CustomPssmInfo" )
		.def_readonly( "name", &custom_pssm::name )
		.def_readonly( "counts", &custom_pssm::counts )
		;

	register_ptr_to_python< custom_pssm::ptr >();

	def( "get_all_custom_pssms", get_all_custom_pssms );
	def( "get_custom_pssms", get_custom_pssms );
	def( "get_custom_pssm_set_names", get_custom_pssm_set_names );
	def( "get_custom_pssm", get_custom_pssm );

}



} //namespace biopsy
