/**
 * Copyright John Reid 2013
 *
 * @file Code to test BiFa analysis maximal chain.
 */

#define BOOST_TEST_MODULE bifa_score
#include <boost/test/unit_test.hpp>

#include <biopsy/init.h>
#include <biopsy/analyse.h>
#include <biopsy/pssm.h>

BOOST_AUTO_TEST_CASE( test_bifa_score )
{
    using namespace biopsy;

    init();
    pssm_parameters::singleton().use_score = false;

    string_vec_ptr pssm_names( new string_vec );
    pssm_names->push_back( "M00023" );
    pssm_names->push_back( "M00436" );
    pssm_names->push_back( "M00716" );
    pssm_names->push_back( "M00938" );
    pssm_names->push_back( "R04653" );

    sequence seq = "ACGCGAGCAGGGTCATTAAATCGAGCGTCGCGGCGCGCGACAAGGACGGCATTATTAGCGTGCTACGACTACGACTTG";

    binding_hit::vec_ptr hits = score_pssms_on_sequence( pssm_names, seq );

	namespace fs = boost::filesystem;
	biopsy::build_svg(
		fs::path( "test-bifa.svg" ),
		"Test BiFa",
		seq,
		0.03,
		*hits,
		10,
		false,
		false,
		0,
		"",
		1.
	);
}
