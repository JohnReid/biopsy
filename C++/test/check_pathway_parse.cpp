#include <bio/amigo_pathways.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <string>
#include <iostream>
using namespace std;


void check_pathway_parse()
{
	cout << "******* check_pathway_parse()" << endl;

#ifdef VERBOSE_CHECKING
	progress_timer timer;
#endif


	pathways.parse_default_xml_files();
}

void register_pathway_parse_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_pathway_parse), 0);
}


