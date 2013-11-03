#include "bio_test_data.h"

#include <bio/ott_match_normaliser.h>
#include <bio/matrix.h>
USING_BIO_NS;


#include <boost/progress.hpp>
#include <boost/test/unit_test.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <fstream>
#include <vector>
using namespace std;

#ifdef _DEBUG
# define TEST_SEQ_LENGTH 1000
#else
# define TEST_SEQ_LENGTH 100000
#endif



//#define VERBOSE_CHECKING

void check_ott_normalisation()
{
	ensure_biobase_params_built();

	cout << "******* check_ott_normalisation()" << endl;

	//number of quanta in normalisation
	const size_t num_quanta = 10;

	//create a random sequence
	const size_t test_seq_length = TEST_SEQ_LENGTH;
	seq_t test_seq;
	test_seq.reserve(test_seq_length);
	generate_random_nucleotide_seq(inserter(test_seq, test_seq.begin()), test_seq_length);

	//create a vector with 2 biobase pssms in it
	typedef BiobasePssm<PssmEntry> pssm_t;
	vector <pssm_t> pssms;
	vector <string> names;
	pssm_t pssm1(acgt_pssm->begin(), acgt_pssm->end());
	pssms.push_back(pssm1);
	names.push_back("acgt");
	pssms.push_back(BiobasePssm<PssmEntry>(V$DEAF1_01_pssm->begin(), V$DEAF1_01_pssm->end()));
	names.push_back("V$DEAF1_01");

	//for each pssm normalise the match scores
	for (int i = 0; (size_t) i < pssms.size(); ++i)
	{
#ifdef VERBOSE_CHECKING
		progress_timer timer;
#endif

		OttNormaliser::ptr_t ott_normaliser(new OttNormaliser(pssms[i], test_seq.begin(), test_seq.end(), num_quanta));

#ifdef VERBOSE_CHECKING
		copy(
			ott_normaliser->normalised_scores.begin(),
			ott_normaliser->normalised_scores.end(),
			ostream_iterator<float_t>(cout, ","));
		cout << endl;
#endif

		OttNormaliserCache::singleton().put_normaliser(names[i], ott_normaliser);
	}
}

BOOST_TEST_DONT_PRINT_LOG_VALUE( OttNormaliserCache )


void check_ott_normaliser_cache() {

	cout << "******* check_ott_normaliser_cache()" << endl;

	// create and open a character archive for output
    std::ofstream ofs("ott_normaliser_cache.txt");
    boost::archive::text_oarchive oa(ofs);

    // write class instance to archive
    oa << OttNormaliserCache::singleton();
    // close archive
    ofs.close();

    // ... some time later restore the class instance to its orginal state
    // create and open an archive for input
    std::ifstream ifs("ott_normaliser_cache.txt", std::ios::binary);
    boost::archive::text_iarchive ia(ifs);
    // read class state from archive
    OttNormaliserCache new_cache;
    ia >> new_cache;
    // close archive
    ifs.close();

	BOOST_CHECK_EQUAL(new_cache, OttNormaliserCache::singleton());
}


