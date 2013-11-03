#include <bio/site_data.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <string>
#include <iostream>
using namespace std;

//#define VERBOSE_CHECKING



void check_site_data_parse()
{
	namespace fs = boost::filesystem;

	cout << "******* check_site_data_parse()" << endl;

	SiteData::vector_t site_data;
	fs::ifstream stream(fs::path("output.txt"));
	bool parsed = SiteData::parse(stream, site_data);
}

void
register_site_data_tests(boost::unit_test::test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_site_data_parse), 0);
}
