

#include "bio_test_data.h"

#include <bio/biobase_match.h>
#include <bio/matrix.h>
#include <bio/score_sequence.h>
USING_BIO_NS;

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/progress.hpp>
#include <boost/test/parameterized_test.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <string>
using namespace std;

/**
Results from TRANSFAC for the test sequence:

matrix                    position  core   matrix sequence (always the               factor name
identifier                (strand)  match  match  (+)-strand is shown)

V$DEAF1_01                  86 (+)  1.000  0.857  gcgcctcagtgtacTTCCGaacgaa          DEAF1
V$IPF1_Q4_01               111 (+)  1.000  0.975  tgagtCATTAataga                    IPF1
V$EFC_Q6                   336 (-)  0.910  0.853  ttacccgaGTGACt                     RFX1
                                                                                     (EF-C)
V$AP1_Q4                   343 (+)  1.000  0.990  agTGACTgact                        AP-1
V$CREBP1_01                384 (+)  1.000  1.000  TTACGtaa                           CRE-BP1
V$CREBP1_01                384 (-)  1.000  1.000  ttaCGTAA                           CRE-BP1
V$HNF3ALPHA_Q6             446 (+)  1.000  0.973  TGTTTgctctc                        HNF-3alpha
V$FOXP3_Q4                 589 (+)  0.859  0.831  tatgaGCTGTttcagat                  FOXP3
V$HNF1_Q6                  694 (-)  0.952  0.873  cttctgaaggAGTAActc                 HNF-1
V$COMP1_01                 888 (+)  1.000  0.859  taggaaGATTGgccacatcagctt           COMP1
V$CDPCR1_01                979 (-)  0.896  0.914  acgaTCTATg                         CDP CR1
*/

void check_biobase_matches(BiobaseParams const & params)
{
	ensure_biobase_params_built();

	cout << "******* check_biobase_matches(): Testing with " << params.name.c_str() << " parameters" << endl;
	//progress_timer timer;

	boost::io::ios_precision_saver saver(cout);
	cout.precision(3);

	//build the pssms
	const Pssm matrix_pssm = *params.pssm;
	const Pssm iupac_pssm = make_pssm_from_iupac(params.iupac->begin(), params.iupac->end());

#ifdef VERBOSE_CHECKING
	cout
		<< "Matrix pssm:" << endl
		<< matrix_pssm
		<< endl;
	cout
		<< "IUPAC pssm:" << endl
		<< iupac_pssm
		<< endl;
#endif

	//match to compute the score
	Hit pssm_match (-1.0, 0);
	Hit iupac_match (-1.0, 0);
	if (params.run_over_range)
	{
		pssm_match =
			score_sequence(
				matrix_pssm, 
				params.seq_begin, 
				params.seq_end, 
				params.match_complement);
		iupac_match =
			score_sequence(
				iupac_pssm, 
				params.seq_begin, 
				params.seq_end,
				params.match_complement);
	}
	else
	{
		assert(params.seq_end - params.seq_begin == int( matrix_pssm.size() ) );
		pssm_match.score = matrix_pssm.score(params.seq_begin, params.match_complement);
		iupac_match.score = iupac_pssm.score(params.seq_begin, params.match_complement);
	}

	//compare against expected
    BOOST_CHECK_CLOSE(params.target_match.score, pssm_match.score, 1.0); //only compare to 1% (we have 3 figure accuracy)
	BOOST_CHECK_CLOSE(iupac_match.score, pssm_match.score, params.iupac_tolerance); //IUPAC should be close to PSSM match
    BOOST_CHECK_EQUAL(params.target_match.position, pssm_match.position); //position of the pssm match
    BOOST_CHECK_EQUAL(params.target_match.complement, pssm_match.complement); //was it on the complementary strand
}


void register_biobase_matches_tests(boost::unit_test::test_suite * test)
{
	ensure_biobase_params_built();

	test->add(BOOST_PARAM_TEST_CASE(&check_biobase_matches, biobase_params.begin(), biobase_params.end()), 0);

	//superceded
	//test->add(BOOST_TEST_CASE(&check_ott_normaliser_cache), 0);
	//test->add(BOOST_TEST_CASE(&check_ott_normalisation), 0);
	//test->add(BOOST_TEST_CASE(&check_weak_pssms), 0);
}

