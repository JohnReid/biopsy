
#include <bio/remo.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include <iostream>
using namespace std;


//#define VERBOSE_CHECKING

void ensure_test_remos_built()
{
	static bool already_done = false;
	if (! already_done) {
		build_test_remos();
		already_done = true;
	}
}


void
check_remos()
{
	cout << "******* check_remos()" << endl;
	ensure_test_remos_built();

#ifdef VERBOSE_CHECKING
	cout << test_remos;
#endif
}


void register_remos_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_remos), 0);
}


