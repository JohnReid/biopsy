

#include <bio/random.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;


void
check_random()
{
	cout << "******* check_random()" << endl;

	//make sure that we can seed the random numbers at will
	seed_default_rng(1234);
	BOOST_CHECK_EQUAL(get_uniform_index(5), 3u);
	BOOST_CHECK_EQUAL(get_uniform_index(5), 4u);
	BOOST_CHECK_EQUAL(get_uniform_index(5), 3u);
	BOOST_CHECK_EQUAL(get_uniform_index(5), 3u);
	BOOST_CHECK_EQUAL(get_uniform_index(5), 4u);

	//make sure our random values are within range - check a few times
	for (size_t i = 1; i < 1000; ++i)
	{
		BIO_NS::float_t rnd_01 = BIO_NS::float_t(get_uniform_01());
		BOOST_CHECK(0.0 <= rnd_01);
		BOOST_CHECK(rnd_01 <= 1.0);
	}

	//do tests over various ranges
	for (size_t max = 1; max < 100; ++max)
	{
		//check we hit every number in the range - this _could_ take a long time
		for (size_t num = 0; max != num; ++num)
		{
			size_t rnd_idx;
			while (num != (rnd_idx = get_uniform_index(max)))
			{
				//check we are in the right range
				BOOST_CHECK(rnd_idx < max);
			}
		}
	}
}

void register_random_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_random), 0);
}


