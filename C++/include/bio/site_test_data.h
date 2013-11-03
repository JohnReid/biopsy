#ifndef BIO_SITE_TEST_DATA
#define BIO_SITE_TEST_DATA

#include "bio/defs.h"
#include "bio/site.h"
#include "bio/site_test_case.h"
#include "bio/test_data.h"


BIO_NS_START

struct SiteTestData
	: BiFaTestData
{
	Site * site;

public:
	SiteTestData(
		const SiteTestCase & test_case);

	virtual ~SiteTestData();

	Test::ptr_t get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm );
};


BIO_NS_END

#endif //BIO_SITE_TEST_DATA
