

#include <bio/svg_match.h>
#include <bio/biobase_match.h>
USING_BIO_NS

#include <boost/progress.hpp>
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;
using namespace boost;

#include <xercesc/util/PlatformUtils.hpp>
XERCES_CPP_NAMESPACE_USE

#include <iostream>
using namespace std;

void check_create_match_svg() {

	cout << "******* check_create_match_svg()" << endl;

#ifdef VERBOSE_CHECKING
	progress_timer timer;
#endif


	XMLPlatformUtils::Initialize();

	const BIO_NS::float_t match_threshold = .05f;
	const std::string title = "Test";
	const seq_t seq = "atcagtcgatcagctagctagtcagc";
	const bool show_titles = true;
	SvgDomDocument doc(seq.size(), match_threshold, 1.0f, title, seq, show_titles);

	dom_print(doc.doc->getDocumentElement(), "test.svg");

	XMLPlatformUtils::Terminate();
}


void register_create_match_svg_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_create_match_svg), 0);
}


