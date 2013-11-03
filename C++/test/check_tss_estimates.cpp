/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#include "bio_test_defs.h"
#include "bio/tss_estimates.h"
#include "bio/environment.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/filesystem/operations.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;



const std::vector< std::string > estimate_lines = boost::assign::list_of
	( "0;ENSMUSG00000022762;ENSMUST00000067602;mus_musculus_core_37_34e;5prime;PROOF;chromosome;16;80308349;1;-231995;CCGCCTCCTCCCGCCGGCTCCGCTGCCGCCGCGGCGGCCGCTGCTGCCGC;STANDARDCLONE;HGI;THC2260281;mouse;human;(-231995|-231739),;" )
	( "1;ENSMUSG00000022762;ENSMUST00000067602;mus_musculus_core_37_34e;5prime;PROOF;chromosome;16;80308349;1;-231995;CCGCCTCCTCCCGCCGGCTCCGCTGCCGCCGCGGCGGCCGCTGCTGCCGC;STANDARDCLONE;HGI;THC2260281;mouse;chimp;(-231995|-231739),;" )
	( "2;ENSMUSG00000022762;ENSMUST00000037785;mus_musculus_core_37_34e;5prime;PROOF;chromosome;16;80308349;1;-231995;CCGCCTCCTCCCGCCGGCTCCGCTGCCGCCGCGGCGGCCGCTGCTGCCGC;STANDARDCLONE;HGI;THC2260281;mouse;human;(-231995|-231739),;" )
	( "3;ENSMUSG00000022762;ENSMUST00000037785;mus_musculus_core_37_34e;5prime;PROOF;chromosome;16;80308349;1;-231995;CCGCCTCCTCCCGCCGGCTCCGCTGCCGCCGCGGCGGCCGCTGCTGCCGC;STANDARDCLONE;MGCHUMAN;gi_31324930_gb_BC052946.1_;mouse;human;(-231995|-231739),;" )
	( "4;ENSMUSG00000022762;ENSMUST00000037785;mus_musculus_core_37_34e;5prime;PROOF;chromosome;16;80308349;1;-231995;CCGCCTCCTCCCGCCGGCTCCGCTGCCGCCGCGGCGGCCGCTGCTGCCGC;STANDARDCLONE;HGI;THC2260281;mouse;chimp;(-231995|-231739),;" )
	( "5;ENSMUSG00000022762;ENSMUST00000037785;mus_musculus_core_37_34e;5prime;PROOF;chromosome;16;80308349;1;-231995;CCGCCTCCTCCCGCCGGCTCCGCTGCCGCCGCGGCGGCCGCTGCTGCCGC;STANDARDCLONE;MGCHUMAN;gi_31324930_gb_BC052946.1_;mouse;chimp;(-231995|-231739),;" )
	;

//#define VERBOSE_CHECKING
void
check_tss_estimate_line_parse( const std::string & line )
{
	cout << "******* check_tss_estimate_line_parse()\n";

	TssEstimate estimate;
	estimate.parse( line );
}

void
check_tss_estimate_file_parse()
{
	cout << "******* check_tss_estimate_file_parse()\n";

	//remove serialised copy if already there
	namespace fs = boost::filesystem;
	const fs::path serialised(
		BioEnvironment::singleton().get_serialised_tss_estimates_file()
	);
	if( fs::exists( serialised ) )
	{
		fs::remove( serialised );
	}

	TssEstimates::singleton();
}


void
register_tss_estimates_tests( boost::unit_test::test_suite * test )
{
	test->add( BOOST_TEST_CASE( &check_tss_estimate_file_parse ), 0);
	test->add( BOOST_PARAM_TEST_CASE( &check_tss_estimate_line_parse, estimate_lines.begin(), estimate_lines.end() ), 0);
}


//
// Are we going to compile this test into its own executable?
//
#ifdef BIO_STANDALONE_TEST
BIO_DEFINE_STANDALONE_TEST( "TSS estimates", register_tss_estimates_tests )
#endif //BIO_STANDALONE_TEST


