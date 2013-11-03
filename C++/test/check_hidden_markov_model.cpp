
#include <bio/environment.h>
#include <bio/hidden_markov_model.h>
#include <bio/hmm_gen_sequence.h>
#include <bio/chromosomes_file_set.h>
#include <bio/hmm_forward_backward.h>
#include <bio/hmm_viterbi.h>
#include <bio/hmm_baum_welch.h>
#include <bio/sequence.h>
#include <bio/fasta.h>
#include <bio/hmm_dna.h>
USING_BIO_NS;

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign.hpp>
#include <boost/progress.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <gsl/gsl_randist.h>

#include <iostream>
#include <list>
#include <fstream>
using namespace std;

#include <math.h>


//#define VERBOSE_CHECKING

#define ENSEMBL_DIR "C:\\data\\ensembl"
#define CHROMO_DIR ENSEMBL_DIR "\\chromosomes"
#define FASTA_FILE ENSEMBL_DIR "\\test_seq_5.fa"
#define LARGE_FASTA_FILE ENSEMBL_DIR "\\Mus_musculus.NCBIM33.may.dna.contig.fa"
//#define LARGE_FASTA_FILE ENSEMBL_DIR "\\test_seq_1.fa"

#define PROB_TOLERANCE 0.01 //give about 1% tolerance on prob comparisons

#define NUM_BAUM_WELCH_ITERATIONS 10


bool
new_prob_less(
	prob_t last_prob,
	prob_t this_prob)
{
	//make sure we only have small positive log probs (due to rounding errors)
	if (0 < last_prob && last_prob < 1.0e-015) {
		last_prob = 0.0;
	}
	return last_prob <= this_prob + PROB_TOLERANCE * abs(this_prob) / 2;
}


typedef MarkovState<NucleoCode> state_t;
typedef HiddenMarkovModel<NucleoCode> hmm_t;
typedef HmmSequenceGenerator<NucleoCode> seq_gen_t;
BOOST_TEST_DONT_PRINT_LOG_VALUE(hmm_t);

hmm_t hmm;

const SeqList test_seqs = boost::assign::list_of
	("AT")
	("ATC")
	("")
	("A")
	("C")
	("G")
	("T")
	("AA")
	("CC")
	("ACGCTCGT")
	("CTTCTTCTCTCTCTAGAGGAGAGAGA")
	("CTTCTAGAGGAGTCTCTCTCTAGAGA")
	("CTTCTTCTCTCTCT")
	("ATCATCATCATCATCATCATCATCATCATCATCATCATC")
	("CTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCT")
	("CTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCT")
	("CTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTCTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCTAGAGGAGTCTCTCTCTAGAGACTTCT")
	("CGTGGATGGGTCCTC")
	("TTCTTA")
	("ATAATTAAA")
	("TAAATTTTA")
	("TGACG")
	("ATTAATTT")
	("ACTTGGTAACA")
	("TATAAATTT")
	("TAGGAAACAAA")
	("AATAT")
	("GTTCAACAGG")
	("CGCAGGCC")
	("TTAATTTATAAAAAA")
	("ATCGCGGA")
	("TAAAATATAA")
	("CAGAATGCTATGAC")
	("AGTT")
	("CGA")
	("AAATTTTAATTATA")
	("GGAGTGCCAA")
	("CCGACAGC")
	("CCATAATGGG")
	("TGCCAGGGGGTCTGAC")
	("GGTTCTTGATAGT")
	("CCAG")
	("GCTTTCCGT")
	("ACAAAGGAAC")
	("TGCACCCC")
	("ATTTAATATTT")
	("TATTCCGCTG")
	("ACTGCAT")
	("AAATTAAT")
	("GGCCGGGAA")
	("CAGCTGATTTCCTA")
	("ATATATAAT")
	("TTATATTTT")
	("ACCGCC")
	("CACACGCAGG")
	("ATCGGACCA")
	("ATATTTTGGC")
	("CCTGGTTAGAGAG")
	("TAGACTCGATCGAA")
	("ACCGGGAGGGCCGCTCC")
	("TGTCATGTCC")
	("ACACCCGAAGAGTA")
	("AGGACATGTC")
	("AAGGGAGGATTA")
	("AGTTGACCCTA")
	("CGAGGTCGGGC")
	("AGGGTGCTGAG")
	("TGCCCCCG")
	("ATTTATTA")
	("CCTCCTGGCT")
	("GCTTGGGGCC")
	("TAGAGGTAC")
	("GTTGCGCATCG")
	("ATAAAAAATTTA")
	("ATGCAGTC")
	("CAAGATACG")
	("TTTTATTTAATTA")
	("CTAAGCGA")
	("CCATGATCGCT")
	("ATTAATTTT")
	("AGAGAGAGCGAC")
	("TTTCCGTCTC")
	("ATAATTATTTA")
	("TAATTTATTAA")
	("TTTTTGTG")
	("AAGATGGTGTT")
	("TTATTATA")
	("AAATTTATATTTA")
	("TTCTGGCTTA")
	("TCAGATTAC")
	("ATTAATAATTTTTTAA")
	("GCCCAGGCCGGGT")
	("AAATAATT")
	("GGCT")
	("AAAAAAAAATTAA")
	("GCCCGCCGG")
	("CCACAAGA")
	("TGAAGGTTATCG")
	("TCTTGGAGAAAA")
	("TCCTGATTCTATAG")
	("TAAAATATAATTAAA")
	("AAAATA")
	("TTAACC")
	("TATAATAAAAT")
	("TTTAAAATT")
	("GCATGGGGCGT")
	("GTAACCTGTGGCG")
	("GCCCGGGG")
	("ATAAAGACGTCCGTC")
	("AATTCCTGAG")
	("AAAGACACCGAAA")
	("GGGGACCGC")
	("TTTTTAAATAA")
	("TATAAAAATT")
	("CATCGCCACGCAAGC")
	("TCGGGGTG")
	("GTCAGCT")
	;


typedef list<size_t> index_list_t;
BOOST_TEST_DONT_PRINT_LOG_VALUE(index_list_t);

seq_t long_test_seq;

typedef vector<SeqList> seq_vec_vec;
const seq_vec_vec hmm_multiple_seqs = boost::assign::list_of< SeqList >
	(boost::assign::list_of< seq_t >
		(seq_t(""))
		(seq_t(""))
		(seq_t(""))
		)
	(boost::assign::list_of< seq_t >
		(seq_t("A"))
		(seq_t("T"))
		(seq_t(""))
		)
	(boost::assign::list_of< seq_t >
		(seq_t(""))
		(seq_t("T"))
		(seq_t("TC"))
		)
	(boost::assign::list_of< seq_t >
		(seq_t("ACGT"))
		(seq_t("TGCA"))
		)
	(boost::assign::list_of< seq_t >
		(seq_t("ATATATATATATAAATATATAATATATATTTTAA"))
		(seq_t("ACACACCCACCACCCACACACAAAACACACACCC"))
		(seq_t("GTGTGTGTTGTGTGTGGTGTTTTTGTGTGGTGTG"))
		)
	(boost::assign::list_of< seq_t >
		(seq_t("ATCGACTAGCAGCGAGT"))
		(seq_t("GATAGATAGAGAGTA"))
		(seq_t("CTGCTGCTCTTCG"))
		)
	;


/** Functor to check the Baum-Welch algorithm for a model on a sequence. */
struct CheckBaumWelch : std::binary_function<hmm_t, seq_t, void>
{
	void operator()(hmm_t hmm, const seq_t & sequence) const
	{
		const prob_t likelihood_before_training =
			get_per_char_likelihood(hmm, sequence.begin(), sequence.end());

		for (unsigned i = 0; NUM_BAUM_WELCH_ITERATIONS != i; ++i)
		{
			baum_welch_single(
				hmm,
				sequence.begin(),
				sequence.end(),
				sequence.rbegin(),
				sequence.rend());
		}

		const prob_t likelihood_after_training =
			get_per_char_likelihood(hmm, sequence.begin(), sequence.end());

#ifdef VERBOSE_CHECKING
		std::cout
			<< "Size: " << hmm.states.size() << "; "
			<< likelihood_before_training << ", "
			<< likelihood_after_training << '\n';
#endif

		//check we are not more than 5% worse than previously
		BOOST_CHECK(likelihood_after_training >= likelihood_before_training - 0.05);
	}
};


/** Functor to check baum welch on a list of sequences. */
struct CheckBaumWelchOnSequences
{
	const SeqList & sequences;
	CheckBaumWelchOnSequences(const SeqList & sequences)
		: sequences(sequences)
	{ }
	void operator()(hmm_t & hmm) const
	{
		std::for_each(sequences.begin(), sequences.end(), std::bind1st(CheckBaumWelch(), hmm));
	}
};



/** Print some random test sequences to cout. */
void
generate_test_sequences()
{
	cout << "******* generate_test_sequences()" << endl;

	const unsigned mean_length = 10;
	const unsigned num_sequences = 100;

	gsl_rng * r;
	/* create a generator chosen by the 
		environment variable GSL_RNG_TYPE */
	{
		const gsl_rng_type * T;
	
		gsl_rng_env_setup();

		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
	}

	typedef boost::array<prob_t, 3> base_dist_t;

	const base_dist_t uniform_dist		= boost::assign::list_of( 0.25 )( 0.50 )( 0.75 );
	const base_dist_t a_t_dist		= boost::assign::list_of( 0.50 )( 0.50 )( 0.50 );
	const base_dist_t c_g_dist		= boost::assign::list_of( 0.10 )( 0.50 )( 0.90 );

	for (unsigned s = 0; num_sequences != s; ++s)
	{
		//choose a dist
		const prob_t dist_rnd = get_uniform_01();
		const base_dist_t & dist = dist_rnd < .5 ? uniform_dist : (dist_rnd < .75 ? a_t_dist : c_g_dist);

		//choose a length
		const unsigned length = gsl_ran_poisson(r, mean_length);

		//generate a sequence
		seq_t seq;
		for (unsigned i = 0; length != i; ++i)
		{
			const prob_t base_rnd = get_uniform_01();
			if (base_rnd < dist[0])
			{
				seq.push_back('A');
			}
			else if (base_rnd < dist[1])
			{
				seq.push_back('C');
			}
			else if (base_rnd < dist[2])
			{
				seq.push_back('G');
			}
			else
			{
				seq.push_back('T');
			}
		}

		//print the sequence
		cout << "\t(\"" << seq << "\")\n";
	}
	gsl_rng_free(r);
}




void
create_hmm() {

	//1st state
	{
		const prob_t emission_probs [] = { 0.9f, 0.025f, 0.025f, 0.05f };
		prob_vector_t trans_probs;
		trans_probs.push_back(0.1f);
		trans_probs.push_back(0.8f);
		trans_probs.push_back(0.1f);
		hmm.states.push_back(
			state_t(
				0.9f,
				emission_probs,
				trans_probs.begin(),
				trans_probs.end()));
	}

	//2nd state
	{
		const prob_t emission_probs [] = { 0.025f, 0.025f, 0.025f, 0.925f };
		prob_vector_t trans_probs;
		trans_probs.push_back(0.1f);
		trans_probs.push_back(0.1f);
		trans_probs.push_back(0.8f);
		hmm.states.push_back(
			state_t(
				0.05f,
				emission_probs,
				trans_probs.begin(),
				trans_probs.end()));
	}

	//3rd state
	{
		const prob_t emission_probs [] = { 0.025f, 0.9f, 0.025f, 0.05f };
		prob_vector_t trans_probs;
		trans_probs.push_back(0.8f);
		trans_probs.push_back(0.1f);
		trans_probs.push_back(0.1f);
		hmm.states.push_back(
			state_t(
				0.05f,
				emission_probs,
				trans_probs.begin(),
				trans_probs.end()));
	}

#ifdef VERBOSE_CHECKING
	cout << hmm;
#endif
}

void
create_hmm_test_seqs()
{
	ifstream file(FASTA_FILE);
	if (! file.good())
	{
		throw std::logic_error(BIO_MAKE_STRING("Could not open FASTA file: \"" << FASTA_FILE << "\""));
	}

#ifdef VERBOSE_CHECKING
	cout << "Parsing fasta file " FASTA_FILE << endl;
	progress_timer timer;
#endif

	stringstream sequence;
	parse_fasta(file, sequence);

	long_test_seq = sequence.str();
}

void ensure_hmm_built() {
	static bool already_done = false;
	if (! already_done) {
		create_hmm();
		already_done = true;
	}
}

void ensure_hmm_test_seqs_created() {
	static bool already_done = false;
	if (! already_done) {
		create_hmm_test_seqs();
		already_done = true;
	}
}






void check_hmm_overfitting()
{
	ensure_hmm_built();

	cout << "******* check_hmm_overfitting()" << endl;

	//first we overfit the hmm on the first sequence then try to train on the second
	seq_t overfitted_seq("ATCATCATCATCATCATCATCATCATCATCATCATCATC");
	seq_t new_seq("ACTACTACTACTACTACTACTACTACTACTACTACTACT");

	BaumWelchAlgorithm<true> baum_welch;

	//use a new HMM
	hmm_t hmm_2(hmm);

	//first train on the over fitting sequence
	{
		baum_welch.forward(
			hmm_2,
			overfitted_seq.begin(),
			overfitted_seq.end());

		prob_t last_prob = baum_welch.get_log_probability();

#ifdef VERBOSE_CHECKING
		cout
			<< "Overfitted sequence has log probability "
			<< last_prob
			<< " before baum welch "
			<< endl;
#endif

		BOOST_CHECK(0 != BIO_FINITE(last_prob));

		const size_t iterations = NUM_BAUM_WELCH_ITERATIONS;
		for (size_t j = 0; j < iterations; ++j)
		{
			baum_welch.run(
				hmm_2,
				overfitted_seq.begin(),
				overfitted_seq.end(),
				overfitted_seq.rbegin(),
				overfitted_seq.rend());

			BOOST_CHECK(hmm_2.is_consistent());
			//cout << hmm_2;

			baum_welch.forward(
				hmm_2,
				overfitted_seq.begin(),
				overfitted_seq.end());
			const prob_t this_prob = baum_welch.get_log_probability();
			BOOST_CHECK(0 != BIO_FINITE(this_prob));

			BOOST_CHECK(new_prob_less(last_prob, this_prob)); 
			if (! new_prob_less(last_prob, this_prob))
			{
				cerr << "Iteration: " << j << "; " << last_prob << " > " << this_prob << endl;
			}

			last_prob = this_prob;
		}
	} //training on overfitting seq

	//now train on new sequence
	{
		const size_t iterations = NUM_BAUM_WELCH_ITERATIONS;
		for (size_t j = 0; j < iterations; ++j) {

			baum_welch.run(
				hmm_2,
				new_seq.begin(),
				new_seq.end(),
				new_seq.rbegin(),
				new_seq.rend());

			BOOST_CHECK(hmm_2.is_consistent());
			//cout << hmm_2;
		}
	}
}






void check_forward_backward() {

	ensure_hmm_built();
	ensure_hmm_test_seqs_created();

	cout << "******* check_forward_backward(): " << test_seqs.size() << " artificial sequences" << endl;

	for (SeqList::const_iterator i = test_seqs.begin(); test_seqs.end() != i; ++i) {

		const size_t num_obs = i->end() - i->begin();

		ForwardBackwardAlgorithm<true> fw_with_scaling;
		fw_with_scaling.forward(
			hmm,
			i->begin(),
			i->end());
		fw_with_scaling.backward(
			hmm,
			i->rbegin(),
			i->rend());

		ForwardBackwardAlgorithm<false> fw_without_scaling;
		fw_without_scaling.forward(
			hmm,
			i->begin(),
			i->end());
		fw_without_scaling.backward(
			hmm,
			i->rbegin(),
			i->rend());

#ifdef VERBOSE_CHECKING
		cout
			<< "Sequence " << (i - test_seqs.begin())
			<< " (length=" << i->size()
			<< ") has probability (with    scaling): "
			<< fw_with_scaling.get_probability()
			<< endl;
		cout
			<< "Sequence " << (i - test_seqs.begin())
			<< " (length=" << i->size()
			<< ") has probability (without scaling): "
			<< fw_without_scaling.get_probability()
			<< endl;

		cout
			<< "Sequence " << (i - test_seqs.begin())
			<< " (length=" << i->size()
			<< ") has log probability (with    scaling): "
			<< fw_with_scaling.get_log_probability()
			<< endl;
		cout
			<< "Sequence " << (i - test_seqs.begin())
			<< " (length=" << i->size()
			<< ") has log probability (without scaling): "
			<< fw_without_scaling.get_log_probability()
			<< endl;
		cout << endl;
#endif

		//the log probability may underflow to log(0.0) without scaling - allow this in the tests
		if (BIO_FINITE(fw_without_scaling.get_log_probability())) {
			//only do comparison if no underflow
			BOOST_CHECK_CLOSE(fw_with_scaling.get_log_probability(), fw_without_scaling.get_log_probability(), 1.0);
		}
		//in the normal probability case both or neither underflow
		BOOST_CHECK_CLOSE(fw_with_scaling.get_probability(), fw_without_scaling.get_probability(), 1.0);
	}
}










void
check_viterbi()
{
	ensure_hmm_built();
	ensure_hmm_test_seqs_created();

	cout << "******* check_viterbi(): " << test_seqs.size() << " artificial sequences" << endl;

	for (SeqList::const_iterator i = test_seqs.begin(); test_seqs.end() != i; ++i) {

		const size_t num_obs = i->end() - i->begin();

		//doesn't work well without logs on long sequences
		if (num_obs > 20)
		{
			continue;
		}

		index_list_t state_indices;
		ViterbiAlgorithm<false> viterbi_without_logs;
		viterbi_without_logs.viterbi(
			hmm,
			num_obs,
			i->begin(),
			i->end(),
			front_inserter(state_indices));

#ifdef VERBOSE_CHECKING
		cout
			<< "Sequence " << (i - test_seqs.begin())
			<< " has most likely state sequence (not using logs): ";
		for (index_list_t::const_iterator j = state_indices.begin(); state_indices.end() != j; ++j) {
			cout << *j;
		}
		cout << endl;
#endif

		index_list_t state_indices_using_logs;
		ViterbiAlgorithm<true> viterbi_with_logs;
		viterbi_with_logs.viterbi(
			hmm,
			num_obs,
			i->begin(),
			i->end(),
			front_inserter(state_indices_using_logs));

#ifdef VERBOSE_CHECKING
		cout
			<< "Sequence " << (i - test_seqs.begin())
			<< " has most likely state sequence (    using logs): ";
		for (index_list_t::const_iterator j = state_indices_using_logs.begin(); state_indices_using_logs.end() != j; ++j) {
			cout << *j;
		}
		cout << endl;
#endif

		BOOST_CHECK_EQUAL(state_indices, state_indices_using_logs);

#ifdef VERBOSE_CHECKING
		cout << endl;
#endif
	}
}









/** Check the Baum-Welch algorithm on the long sequence. */
void
check_long_test_seq()
{
	ensure_hmm_built();
	ensure_hmm_test_seqs_created();

	cout << "******* check_long_test_seq(): Checking baum welch algorithm on homo sapiens sequence of length " << long_test_seq.size() << endl;

	BaumWelchAlgorithm<true> baum_welch;

	//use a new HMM for each sequence
	hmm_t hmm_2(hmm);

	baum_welch.forward(
		hmm_2,
		long_test_seq.begin(),
		long_test_seq.end());
	prob_t last_prob = baum_welch.get_log_probability();

#ifdef VERBOSE_CHECKING
	cout
		<< "Long test sequence has log probability "
		<< last_prob
		<< " before baum welch "
		<< endl;
#endif

	BOOST_CHECK(0 != BIO_FINITE(last_prob));

	const size_t iterations = NUM_BAUM_WELCH_ITERATIONS;
	for (size_t j = 0; j < iterations; ++j) {

		baum_welch.run(
			hmm_2,
			long_test_seq.begin(),
			long_test_seq.end(),
			long_test_seq.rbegin(),
			long_test_seq.rend());

		BOOST_CHECK(hmm_2.is_consistent());
		//cout << hmm_2;

		baum_welch.forward(
			hmm_2,
			long_test_seq.begin(),
			long_test_seq.end());
		const prob_t this_prob = baum_welch.get_log_probability();
#ifdef VERBOSE_CHECKING
		cout << "Likelihood " << this_prob << endl;
#endif
		BOOST_CHECK(0 != BIO_FINITE(this_prob));

		if (! new_prob_less(last_prob, this_prob)) {
			cerr << last_prob << " > " << this_prob << endl;
		}
		BOOST_CHECK(new_prob_less(last_prob, this_prob));

		last_prob = this_prob;
	}

#ifdef VERBOSE_CHECKING
	cout
		<< "Long test sequence has log probability "
		<< last_prob
		<< " after baum welch "
		<< endl;
#endif
}






/** Check the Baum-Welch algorithm. */
void
check_baum_welch_multiple(const SeqList & seqs)
{
	ensure_hmm_built();

	cout << "******* check_baum_welch_multiple()" << endl;

	//generate the combined sequence
	const seq_t combined_seq = std::accumulate(seqs.begin(), seqs.end(), string(""));

	//train a hmm on the combined sequence
	hmm_t hmm_single(hmm);
	baum_welch_single(
		hmm_single,
		combined_seq.begin(),
		combined_seq.end(),
		combined_seq.rbegin(),
		combined_seq.rend());
	const prob_t single_log_prob = 
		ForwardBackwardAlgorithm<true>()
			.forward(
				hmm_single,
				combined_seq.begin(),
				combined_seq.end()).get_log_probability();


	hmm_t hmm_multiple(hmm);
	prob_t last_multiple_log_prob =
		get_multiple_log_likelihood(hmm_multiple, seqs.begin(), seqs.end());

	//train multiply several times
	for (unsigned i = 0; NUM_BAUM_WELCH_ITERATIONS != i; ++i)
	{
		//train a hmm on all sequences concurrently
		baum_welch_multiple(
			hmm_multiple,
			seqs.begin(),
			seqs.end());

		const prob_t new_log_prob =
			get_multiple_log_likelihood(hmm_multiple, seqs.begin(), seqs.end());
	
#ifdef VERBOSE_CHECKING
		if (new_log_prob < last_multiple_log_prob)
		{
			cout << new_log_prob << " < " << last_multiple_log_prob << endl;
		}
#endif

		last_multiple_log_prob = new_log_prob;
	}

#ifdef VERBOSE_CHECKING
	cout
		<< "Log prob of combined seq after training one at a time: "
		<< single_log_prob
		<< endl;
	cout
		<< "Log prob of combined seq after training all at once: "
		<< last_multiple_log_prob
		<< endl;
#endif //VERBOSE_CHECKING

}


/** Check the Baum-Welch algorithm. */
void
check_baum_welch(const seq_t & seq)
{
	cout << "******* check_baum_welch(): \"" << seq << "\"\n";

	vector<hmm_t> hmms = boost::assign::list_of
		(hmm_t(1))
		(hmm_t(2))
		(hmm_t(3))
		(hmm_t(4))
		(hmm_t(5))
		(hmm_t(6))
		;

	std::for_each(hmms.begin(), hmms.end(), std::bind2nd(CheckBaumWelch(), seq));
}


/** Check the number of states improves learning. */
void
check_more_states_improves_learning(const seq_t & seq)
{
	cout << "******* check_more_states_improves_learning(): \"" << seq << "\"\n";

	typedef vector<hmm_t> hmm_vec;
	hmm_vec hmms = boost::assign::list_of
		(hmm_t(1))
		(hmm_t(2))
		(hmm_t(3))
		(hmm_t(4))
		(hmm_t(5))
		(hmm_t(6))
		;

	prob_t last_likelihood = 0;
	for (hmm_vec::iterator h = hmms.begin(); hmms.end() != h; ++h)
	{
		for (unsigned i = 0; 20 != i; ++i)
		{
			baum_welch_single(
				*h,
				seq.begin(),
				seq.end(),
				seq.rbegin(),
				seq.rend());
		}

		const prob_t new_likelihood = get_per_char_likelihood(*h, seq.begin(), seq.end());

		cout << h->states.size() << ": " << new_likelihood << endl;
		BOOST_CHECK(new_likelihood >= 0.99 * last_likelihood);
		last_likelihood = new_likelihood;
	}
}



void register_hidden_markov_model_tests(test_suite * test)
{
	ensure_hmm_built();
	ensure_hmm_test_seqs_created();

	test->add(BOOST_TEST_CASE(&check_forward_backward), 0);
	test->add(BOOST_PARAM_TEST_CASE(&check_baum_welch, test_seqs.begin(), test_seqs.end()), 0);
	test->add(BOOST_TEST_CASE(&check_hmm_overfitting), 0);
	test->add(BOOST_TEST_CASE(&check_long_test_seq), 0);
	test->add(BOOST_PARAM_TEST_CASE(&check_baum_welch_multiple, hmm_multiple_seqs.begin(), hmm_multiple_seqs.end()), 0);

	//it is not clear how effective these tests are
	//test->add(BOOST_PARAM_TEST_CASE(&check_more_states_improves_learning, test_seqs.begin(), test_seqs.end()), 0);
	//test->add(BOOST_TEST_CASE(&check_viterbi), 2); //can get 2 underflow problems on long sequences when not using logs
	//test->add(BOOST_TEST_CASE(&generate_test_sequences), 0);
}

