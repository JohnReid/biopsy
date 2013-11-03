

#include <bio/test_data.h>
#include <bio/biobase_db.h>
#include <bio/matrix_test_data.h>
#include <bio/biobase_data_traits.h>
USING_BIO_NS;


#include <boost/test/unit_test.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/progress.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <fstream>
#include <vector>
using namespace std;

#define VERBOSE_CHECKING


static const double min_threshold = 0.005;


void
run_tests( BiFaTestData::vec_t & data, const std::string & name )
{
#ifdef VERBOSE_CHECKING
	std::cout << "Running " << data.size() << " " << name << " tests\n";
#endif //VERBOSE_CHECKING

	Test::vec_t tests;
	std::transform(
		data.begin(),
		data.end(),
		back_inserter( tests ),
		BiFaTestData2Test( BiFaAlgorithm::get_default_algorithm() ) );

	test_result_map_t results;
	run_tests(
		tests,
		results,
		min_threshold,
		1.0,
		10 );

#ifdef VERBOSE_CHECKING
	std::cout << name << " test results\n" << results;
#endif //VERBOSE_CHECKING
}

void check_run_tests()
{
	cout << "******* check_run_tests()" << endl;

	static const unsigned seq_length = 200;
	BiFaAlgorithm::ptr_t algorithm(
		new TransfacBiFaAlgorithm( 
			PssmMatchArgs( BIO_NS::float_t( min_threshold ) ),
			true ) );

	{
		BiFaTestData::vec_t tests = boost::assign::list_of
			( BiFaTestData::ptr_t(
				new MatrixTestData(
					BiobaseDb::singleton().get_entry< MATRIX_DATA >( 924 ),
					seq_length ) ) )
			;
		run_tests(
			tests,
			"Single matrix" );
	}

	{
		BiFaTestData::vec_t tests;
		create_test_data_from_sites( tests );
		run_tests(
			tests,
			"Site" );
	}

	{
		BiFaTestData::vec_t tests;
		create_test_data_from_fragments( tests );
		run_tests(
			tests,
			"Fragment" );
	}

	{
		BiFaTestData::vec_t tests;
		create_test_data_from_matrices( tests, seq_length );
		run_tests(
			tests,
			"Matrix" );
	}
}


void register_bifa_test_case_tests(test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_run_tests ), 0 );
}


