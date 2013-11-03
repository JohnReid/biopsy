
#include <bio/sequence.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
#include <bio/biobase_parse_spirit.h>
USING_BIO_NS;
using namespace BIO_NS::spirit;

#include <boost/progress.hpp>
#include <boost/test/unit_test.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <fstream>
using namespace std;


//#define VERBOSE_CHECKING



BOOST_TEST_DONT_PRINT_LOG_VALUE( Matrix::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Site::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Factor::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Fragment::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Molecule::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Pathway::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Compel::map_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( Evidence::map_t );

template <TransData type>
void
check_biobase_table_parse()
{
	typedef DataTraits<type> traits_t;

	cout << "******* check_parsing(): " << traits_t::get_name() << endl;

	//parse the biobase file
	typename traits_t::entry_t::map_t map;
	{
#ifdef VERBOSE_CHECKING
		cout << "Parsing " << traits_t::get_name() << " entries from " << traits_t::get_biobase_file() << endl;
		progress_timer timer;
#endif

		parse<type>(map);

#ifdef VERBOSE_CHECKING
		cout << "Parsed " << map.size() << " entries" << endl;
#endif

	}
	BOOST_CHECK_EQUAL(map.size(), traits_t::get_num_data());
}

//make sure the template functions are instantiated
void check_parse_all_tables()
{
	check_biobase_table_parse< SITE_DATA >();
	check_biobase_table_parse< MATRIX_DATA >();
	check_biobase_table_parse< FACTOR_DATA >();
	check_biobase_table_parse< FRAGMENT_DATA >();
	check_biobase_table_parse< COMPEL_DATA >();
	check_biobase_table_parse< EVIDENCE_DATA >();
	check_biobase_table_parse< PATHWAY_DATA >();
	check_biobase_table_parse< MOLECULE_DATA >();
}



void check_biobase_load_all()
{
	cout << "******* check_biobase_load_all()\n";

	BiobaseDb::singleton().load_all();
}



void
register_biobase_parse_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< MATRIX_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< FACTOR_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< SITE_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< FRAGMENT_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< COMPEL_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< EVIDENCE_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< PATHWAY_DATA > ), 0);
	test->add( BOOST_TEST_CASE( &check_biobase_table_parse< MOLECULE_DATA > ), 0);

	test->add( BOOST_TEST_CASE( &check_biobase_load_all ), 0);

#if 0
#endif
}
