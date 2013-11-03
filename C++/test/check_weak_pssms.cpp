#include <bio/environment.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
#include <bio/biobase_match.h>
#include <bio/score_sequence.h>
#include <bio/ott_match_normaliser.h>
#include <bio/remo.h>
USING_BIO_NS;

#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

//#define VERBOSE_CHECKING


void
check_weak_pssms()
{
	build_test_remos();

	cout << "******* check_weak_pssms()" << endl;

	Site * site = BiobaseDb::singleton().get_entry<SITE_DATA>(TableLink(SITE_DATA, 2156));
	BOOST_CHECK_EQUAL(site->get_name(), std::string("GATA1$CONS_01")); //if database is updated this could change

	BiobasePssm<IupacCode> pssm(site->sequence.begin(), site->sequence.end());

	const size_t norm_seq_length = 10000;
	seq_t norm_seq;
	norm_seq.reserve(norm_seq_length);
	generate_random_nucleotide_seq(inserter(norm_seq, norm_seq.begin()), norm_seq_length);
	OttNormaliser normaliser(pssm, norm_seq.begin(), norm_seq.end());

	//match the pssm against a remo
	typedef std::vector<Hit> result_vec_t;
	result_vec_t results;
	const seq_t & remo_seq = test_remos.begin()->second->map[HUMAN_SPECIES];
	score_sequence(
		pssm,
		remo_seq.begin(),
		remo_seq.end(),
		true,
		inserter(results, results.end()));

#ifdef VERBOSE_CHECKING
	copy(results.begin(), results.end(), ostream_iterator<Hit>(cout, "\n"));
	cout << endl << endl;
#endif

	normaliser.normalise(results.begin(), results.end());

#ifdef VERBOSE_CHECKING
	copy(results.begin(), results.end(), ostream_iterator<Hit>(cout, "\n"));
	cout << endl << endl;
#endif
}
