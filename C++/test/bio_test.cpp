/**
@file

Copyright John Reid 2006, 2007, 2013

Contains code to test various functions in the bio library.
*/

#include <bio/options.h>
#include <bio/exceptions.h>
#include <bio/gsl.h>
USING_BIO_NS;


#include <antlr/ANTLRException.hpp>
using namespace antlr;


#include <boost/progress.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_monitor.hpp>
#include <boost/test/parameterized_test.hpp>
using namespace boost;
using namespace boost::unit_test;


#include <xercesc/dom/DOMException.hpp>
XERCES_CPP_NAMESPACE_USE



#include <iostream>
#include <iterator>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;


void register_adjust_phylo_tests( test_suite * test );
void register_bifa_test_case_tests( test_suite * test );
void register_binding_model_tests( test_suite * test );
void register_binding_max_chain_tests( test_suite * test );
void register_biobase_likelihoods_tests( test_suite * test );
void register_biobase_matches_tests( test_suite * test );
void register_biobase_parse_tests( test_suite * test );
void register_cache_tests( test_suite * test );
void register_clustering_tests( test_suite * test );
void register_contingency_homogeneity_tests( test_suite * test );
void register_create_match_svg_tests( test_suite * test );
void register_factor_pathway_map_tests( test_suite * test );
void register_factor_synonym_tests( test_suite * test );
void register_fasta_tests( test_suite * test );
void register_hidden_markov_model_tests( test_suite * test );
void register_hmm_dna_tests( test_suite * test );
void register_kegg_tests( test_suite * test );
void register_markov_model_tests( test_suite * test );
void register_math_tests( test_suite * test );
void register_matrix_dependencies_tests( test_suite * test );
void register_matrix_match_map_tests( test_suite * test );
void register_multi_seq_match_tests( test_suite * test );
void register_pathway_parse_tests( test_suite * test );
void register_pssm_likelihoods_tests( test_suite * test );
void register_pssm_match_tests( test_suite * test );
void register_pssm_motif_tests( test_suite * test );
void register_pssm_pathway_map_tests( test_suite * test );
void register_random_tests( test_suite * test );
void register_remo_parse_tests( test_suite * test );
void register_remos_tests( test_suite * test );
void register_score_map_tests( test_suite * test );
void register_site_data_tests( test_suite * test );
void register_site_data_tests( test_suite * test );
void register_serialisation_strategy_tests( test_suite * test );
void register_species_file_sets_tests( test_suite * test );
void register_svg_tests( test_suite * test );
void register_tss_estimates_tests( test_suite * test );
void register_wsdl_tests( test_suite * test );



test_suite*
init_unit_test_suite( int argc, char * argv [] )
{
    test_suite* test = BOOST_TEST_SUITE("Bio test suite");

    try
    {
        BioOptions::singleton().parse(argc, argv);

        gsl_init();

        //ensure_hmm_multiple_seqs_built();

        //register exception handlers
        //unit_test_monitor.register_exception_translator<DOMException>(&translate_dom_exception);
        //unit_test_monitor.register_exception_translator<ANTLRException>(&translate_antlr_exception);
        //unit_test_monitor.register_exception_translator<string>(&translate_string_exception);
        //unit_test_monitor.register_exception_translator<const char *>(&translate_char_exception);

        register_adjust_phylo_tests( test );
        //register_binding_max_chain_tests( test );
        register_binding_model_tests( test );
        register_pssm_match_tests( test );
        //register_kegg_tests( test );
        register_bifa_test_case_tests( test );
        register_tss_estimates_tests( test );
        //register_wsdl_tests( test );
        //register_serialisation_strategy_tests( test );
        register_pssm_pathway_map_tests( test );
        //register_biobase_parse_tests( test );
        register_cache_tests( test );
        register_pssm_likelihoods_tests( test );
        register_biobase_likelihoods_tests( test );
        register_factor_synonym_tests( test );
        register_score_map_tests( test );
        //register_hmm_dna_tests( test );
        register_hidden_markov_model_tests( test );
        register_create_match_svg_tests( test );
        register_matrix_dependencies_tests( test );
        register_pathway_parse_tests( test );
        register_pssm_motif_tests( test );
        register_random_tests( test );
        register_species_file_sets_tests( test );
        register_svg_tests( test );
        register_remos_tests( test );
        register_site_data_tests( test );
        register_matrix_match_map_tests( test );
        register_multi_seq_match_tests( test );
        //register_markov_model_tests( test );
        register_math_tests( test );
        register_biobase_matches_tests( test );
        //register_clustering_tests( test );
        register_contingency_homogeneity_tests( test );
        register_factor_pathway_map_tests( test );
        register_fasta_tests( test );
        register_remo_parse_tests( test );

    }
    catch (const std::exception & e)
    {
        cerr << "Exception: " << e.what() << endl;
    }

    return test;
}

