
#include "bio_test_data.h"

#include <bio/multi_seq_match.h>
#include <bio/remo.h>
USING_BIO_NS;

#include <boost/progress.hpp>
#include <boost/test/unit_test.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;

//#define VERBOSE_CHECKING

typedef std::pair<seq_t::const_iterator, seq_t::const_iterator> seq_def_t;
typedef std::vector<seq_def_t> seqs_vec;

void add_seq(seqs_vec & seqs, const seq_t & seq)
{
	seqs.push_back(make_pair(seq.begin(), seq.end()));
}

void add_seqs(seqs_vec & seqs, const char * remo_name)
{
	add_seq(seqs, test_remos[remo_name]->map[HUMAN_SPECIES]);
	add_seq(seqs, test_remos[remo_name]->map[MOUSE_SPECIES]);
	add_seq(seqs, test_remos[remo_name]->map[RAT_SPECIES]);
	add_seq(seqs, test_remos[remo_name]->map[DOG_SPECIES]);
	add_seq(seqs, test_remos[remo_name]->map[XENOPUS_SPECIES]);
	add_seq(seqs, test_remos[remo_name]->map[TETRAODON_SPECIES]);
}

void check_multiple_sequence_match()
{
	cout << "******* check_multiple_sequence_match()" << endl;

	build_test_remos();

	{
		seqs_vec seqs;
		add_seqs(seqs, "dlx5_1");

		multiple_sequence_match(
			seqs.begin(),
			seqs.end(),
			BIO_NS::float_t(0.999),
			BAYESIAN_SCORE_ALGORITHM);
	}

	{
		seqs_vec seqs;
		add_seqs(seqs, "dlx5_2");

		multiple_sequence_match(
			seqs.begin(),
			seqs.end(),
			BIO_NS::float_t(0.999),
			BAYESIAN_SCORE_ALGORITHM);
	}
}


void register_multi_seq_match_tests(test_suite * test)
{
	ensure_biobase_params_built();

	test->add(BOOST_TEST_CASE(&check_multiple_sequence_match), 0);
}


