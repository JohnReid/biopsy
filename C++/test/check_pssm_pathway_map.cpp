

#include <bio/pathway_associations.h>
USING_BIO_NS;

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include <iostream>
using namespace std;


void
check_pssm_pathway_map()
{
	cout << "******* check_pssm_pathway_map()" << endl;

	BOOST_CHECK_EQUAL( 4252u, PathwayAssociations::singleton().associations.size() );
}


void register_pssm_pathway_map_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_pssm_pathway_map), 0);
}


