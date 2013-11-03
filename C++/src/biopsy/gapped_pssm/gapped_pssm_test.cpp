

#include "biopsy/gapped_pssm.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;


//#define VERBOSE_CHECKING


typedef std::vector< std::string > string_vec;
typedef std::vector< string_vec > string_vec_vec;

string_vec_vec test_seqs = list_of< string_vec >
	(
		list_of< std::string >
			( "gaaaac" )
			( "aacaat" )
			( "aacaat" )
			( "aacaag" )
			( "aaaatc" )
			( "anaacg" )
	)
	(
		list_of< std::string >
			( "ttttt" )
			( "aaaaa" )
	)
	;

void
check_gapped_pssm( const string_vec & test_seqs )
{
	cout << "******* check_gapped_pssm()" << endl;

	using namespace biopsy;
	using namespace biopsy::gapped_pssm;

	dna_vec_list seqs( test_seqs.size() );
	for( unsigned i = 0; test_seqs.size() != i; ++i )
	{
		string_to_dna_vec( test_seqs[i], seqs[i] );
	}

	variational_model model(
		4,
		seqs,
		std::vector< double >( 2, 1.0 ),
		std::vector< double >( 4, 10.0 ),
		std::vector< double >( 4, 0.1 ) );

	double ll = model.log_likelihood();
	for( unsigned i = 0; 25 != i; ++i )
	{
		model.update();
		ll = model.log_likelihood();
	}
	model.update();
}



void
register_compressed_int_array_tests( boost::unit_test::test_suite * test )
{
	test->add( 
		BOOST_PARAM_TEST_CASE( 
			&check_gapped_pssm,
			test_seqs.begin(),
			test_seqs.end() ), 
		0);
}

test_suite*
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite* suite = BOOST_TEST_SUITE( "Bio test suite" );

	try
	{
		register_compressed_int_array_tests( suite );
	}
	catch (const std::exception & e)
	{
		cerr << "Exception: " << e.what() << endl;
	}

    return suite;
}

