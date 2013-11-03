#ifndef BIO_FRAGMENT_TEST_DATA
#define BIO_FRAGMENT_TEST_DATA

#include "bio/defs.h"
#include "bio/fragment.h"
#include "bio/transcription_factor.h"
#include "bio/test_data.h"


BIO_NS_START

struct FragmentTestData
	: BiFaTestData
{
protected:
	const Fragment * fragment;

protected:
	const TF::hit_set_t & get_hits();

public:
	FragmentTestData(
		const Fragment * fragment);

	virtual ~FragmentTestData();

	Test::ptr_t get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm );

};


BIO_NS_END

#endif //BIO_FRAGMENT_TEST_DATA
