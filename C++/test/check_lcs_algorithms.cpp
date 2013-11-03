/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#include "bio_test_defs.h"

#include <bio/lcs_delta_algorithm.h>
#include <bio/remo.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <iostream>
#include <iterator>
#include <list>
using namespace std;

struct WildcardNsEqualTo
{
	bool operator()(char c1, char c2) const
	{
		//'n's are wild cards
		return 'n' == c1 || 'n' == c2 || c1 == c2;
	}
};

void
test_lcs_delta(
	const std::string & s1,
	const std::string & s2,
	const std::string & target)
{
#ifdef VERBOSE_CHECKING
	cout << "Checking LCS(\"" << s1 << "\", \"" << s2 << "\") = \"" << target << "\"" << endl;
#endif

	string result;
	lcs_delta(
		s1.rbegin(),
		s1.rend(),
		s2.rbegin(),
		s2.rend(),
		make_lcs_inserter(result));
	BOOST_CHECK_EQUAL(target, result);
}

void
check_lcs_delta_algorithm()
{
	cout << "******* check_lcs_delta_algorithm()" << endl;

	const std::string nano = "nano";
	const std::string aano = "aano";
	const std::string ano = "ano";
	const std::string anoo = "anoo";
	const std::string nematode_knowledge = "nematode knowledge";
	const std::string empty_bottle = "empty bottle";

	test_lcs_delta(nano, aano, "ano");
	test_lcs_delta(nano, ano, "ano");
	test_lcs_delta(nano, anoo, "ano");
	test_lcs_delta(ano, ano, "ano");
	test_lcs_delta(nano, nematode_knowledge, "nano");
	test_lcs_delta(empty_bottle, nematode_knowledge, "emt ole");

	//use an equals with a wildcard
	{
		const std::string x = "nano";
		const std::string y = "aaao";

		std::string result_x;
		std::string result_y;
		lcs_delta(
			x.rbegin(),
			x.rend(),
			y.rbegin(),
			y.rend(),
			make_lcs_inserter(result_x, result_y),
			WildcardNsEqualTo());
		BOOST_CHECK_EQUAL(y, result_y);
		BOOST_CHECK_EQUAL(x, result_x);
	}

	//test on real sequences
	{
#ifdef VERBOSE_CHECKING
		cout << "Testing DLX5 remo 1" << endl;
#endif

		build_test_remos();

#ifdef VERBOSE_CHECKING
		boost::progress_timer timer;
#endif

		seq_t result;
		lcs_delta(
			test_remos["dlx5_1"]->map[MOUSE_SPECIES].rbegin(),
			test_remos["dlx5_1"]->map[MOUSE_SPECIES].rend(),
			test_remos["dlx5_1"]->map[HUMAN_SPECIES].rbegin(),
			test_remos["dlx5_1"]->map[HUMAN_SPECIES].rend(),
			make_lcs_inserter(result));
		BOOST_CHECK_EQUAL(result.size(), 574u);
	}
}


void register_lcs_algorithms_tests(test_suite * test)
{
	//test->add(BOOST_TEST_CASE(&check_multiple_sequence_match), 0);
	test->add(BOOST_TEST_CASE(&check_lcs_delta_algorithm), 0);
}


//
// Are we going to compile this test into its own executable?
//
#ifdef BIO_STANDALONE_TEST
BIO_DEFINE_STANDALONE_TEST( "lcs_algorithms", register_lcs_algorithms_tests )
#endif //BIO_STANDALONE_TEST

