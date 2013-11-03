
#include <bio/hmm_dna.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/assign.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <sstream>
using namespace std;


//#define VERBOSE_CHECKING


BIO_NS_START

/**
Functor to check conversion.
*/
struct ConversionChecker
{
	template <unsigned order>
	void check(const seq_t & seq) const
	{
		//convert to emissions
		typename DnaHmm<order>::emission_seq_t emissions;
		DnaHmm<order>::convert_to_emission(seq, emissions);

		//convert back again
		seq_t seq_copy;
		DnaHmm<order>::convert_to_dna(emissions, seq_copy);

		//length that should be equal
		unsigned len = (seq.size() / (order + 1)) * (order + 1);
		BOOST_CHECK_EQUAL(emissions.size(), seq.size() / (order + 1));
		BOOST_CHECK_EQUAL(seq_copy.size(), len);
		BOOST_CHECK_EQUAL(seq_copy, seq.substr(0, len));
	}

	void operator()(const seq_t & seq) const
	{
		check<0>(seq);
		check<1>(seq);
		check<2>(seq);
		check<3>(seq);
		check<4>(seq);
		check<5>(seq);
		check<6>(seq);
		check<7>(seq);
		check<8>(seq);
		check<9>(seq);
		check<10>(seq);
		check<11>(seq);
		check<12>(seq);
		check<13>(seq);
		check<14>(seq);
		check<15>(seq);
	}
};

BIO_NS_END



void
check_random_dna_creation()
{
	cout << "******* check_random_dna_creation()" << endl;

	DnaHmmOrderNumStateMap::singleton().insert_model(1, 0, create_dna_model(1, 0));
	DnaHmmOrderNumStateMap::singleton().get_model(1, 0);
	for (unsigned i = 0; 100 != i; ++i)
	{
		seq_t seq;
		DnaHmmOrderNumStateMap::singleton().gen_sequence_from_random_hmm(seq, i);
		BOOST_CHECK_EQUAL(seq.size(), i);

#ifdef VERBOSE_CHECKING
		std::cout << seq << std::endl;
#endif
	}
};

void
check_unsigned_power()
{
	cout << "******* check_unsigned_power()" << endl;

	BOOST_CHECK_EQUAL(1u, unsigned_power(1, 0));
	BOOST_CHECK_EQUAL(1u, unsigned_power(2, 0));
	BOOST_CHECK_EQUAL(1u, unsigned_power(3, 0));
	BOOST_CHECK_EQUAL(1u, unsigned_power(4, 0));

	BOOST_CHECK_EQUAL(0u, unsigned_power(0, 1));
	BOOST_CHECK_EQUAL(1u, unsigned_power(1, 1));
	BOOST_CHECK_EQUAL(2u, unsigned_power(2, 1));
	BOOST_CHECK_EQUAL(3u, unsigned_power(3, 1));
	BOOST_CHECK_EQUAL(4u, unsigned_power(4, 1));

	BOOST_CHECK_EQUAL(0u, unsigned_power(0, 2));
	BOOST_CHECK_EQUAL(1u, unsigned_power(1, 2));
	BOOST_CHECK_EQUAL(4u, unsigned_power(2, 2));
	BOOST_CHECK_EQUAL(9u, unsigned_power(3, 2));
	BOOST_CHECK_EQUAL(16u, unsigned_power(4, 2));

	BOOST_CHECK_EQUAL(0u, unsigned_power(0, 3));
	BOOST_CHECK_EQUAL(1u, unsigned_power(1, 3));
	BOOST_CHECK_EQUAL(8u, unsigned_power(2, 3));
	BOOST_CHECK_EQUAL(27u, unsigned_power(3, 3));
	BOOST_CHECK_EQUAL(64u, unsigned_power(4, 3));
};

void
check_emission_dna_conversion()
{
	cout << "******* check_emission_dna_conversion()" << endl;

	using namespace boost::assign;
	const SeqList test_seqs =
		list_of
			("")
			("A")
			("C")
			("G")
			("T")
			("AAAA")
			("CCCC")
			("GGGG")
			("TTTT")
			("ACGTACGTACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
			("ACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGATACGTACGTACGTATCGATCGAT")
			;

	std::for_each(test_seqs.begin(), test_seqs.end(), ConversionChecker());
}

void
check_hmm_map()
{
	cout << "******* check_hmm_map()" << endl;

	using namespace boost::assign;
	const SeqList training_data =
		list_of
			("acgtacgtacgtatcgatcgat")
			("acgtacgtacgtatcgatcgat")
			("acgtacgtacgtatcgatcgat")
			("acgtacgtacgtatcgatcgat")
			("acgtacgtacgtatcgatcgat")
			;

	DnaHmmOrderNumStateMap hmm_map;
	hmm_map.insert_model(1, 0, create_dna_model(1, 0));
	DnaModel & model = hmm_map.get_model(1, 0);
	model.train(SequenceCollectionList(training_data));
	const prob_t avg_prob_1 = model.get_likelihood(SequenceCollectionList(training_data));
	model.train(SequenceCollectionList(training_data));
	const prob_t avg_prob_2 = model.get_likelihood(SequenceCollectionList(training_data));
	BOOST_CHECK(avg_prob_1 >= avg_prob_2);
}

/** Functor to check training increases sequences' likelihoods. */
struct CheckTrainingIncreasesLikelihoods
{
	typedef std::list<DnaHmmOrderNumStateMap::Index> IndexList;

	unsigned training_iterations;
	const IndexList & hmm_indices;

	CheckTrainingIncreasesLikelihoods(
		const IndexList & hmm_indices,
		unsigned training_iterations = 4)
		: training_iterations(training_iterations)
		, hmm_indices(hmm_indices)
	{
	}

	void operator()(const seq_t & test_seq) const
	{
#ifdef VERBOSE_CHECKING
		std::cout << "\"" << test_seq << "\"\n";
#endif

		DnaHmmOrderNumStateMap hmm_map;

		//train each model
		for (IndexList::const_iterator idx = hmm_indices.begin();
			hmm_indices.end() != idx;
			++idx)
		{

			hmm_map.insert_model(idx->num_states, idx->order, create_dna_model(idx->num_states, idx->order));
			DnaModel & model = hmm_map.get_model(idx->num_states, idx->order);

			const prob_t likelihood_before_training = model.get_likelihood(test_seq);

			for (unsigned i = 0; training_iterations != i; ++i)
			{
				model.train(test_seq);
			}

			const prob_t likelihood_after_training = model.get_likelihood(test_seq);

#ifdef VERBOSE_CHECKING
			std::cout
				<< '(' << idx->num_states << ',' << idx->order << "): "
				<< likelihood_before_training << ", "
				<< likelihood_after_training << '\n';
#endif

			//check we are not more than 5% worse than previously
			BOOST_CHECK(likelihood_after_training >= likelihood_before_training - 0.05);
		}
	}
};

void
check_hmm_single_training()
{
	cout << "******* check_hmm_single_training()" << endl;

	const SeqList test_seqs = boost::assign::list_of
		("AGCTATCAGTCGATGATCGATGCTAGTCG")
		("")
		("c")
		("atatatatatatatatatatatatatatat")
		("AGCTATCAGTCGATGATCGATGCTAGTCG")
		("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
		("cgcgcgcgcgcgcgcgcgcgcgc")
		("agacgattagttagatggcatcgagctaatatcagcagctaatagcgc")
		;

	CheckTrainingIncreasesLikelihoods::IndexList hmm_indices = boost::assign::list_of
		(DnaHmmOrderNumStateMap::Index(1, 0))
		(DnaHmmOrderNumStateMap::Index(2, 0))
		(DnaHmmOrderNumStateMap::Index(3, 0))
		(DnaHmmOrderNumStateMap::Index(1, 1))
		(DnaHmmOrderNumStateMap::Index(2, 1))
		(DnaHmmOrderNumStateMap::Index(3, 1))
		(DnaHmmOrderNumStateMap::Index(1, 2))
		(DnaHmmOrderNumStateMap::Index(2, 2))
		(DnaHmmOrderNumStateMap::Index(3, 2))
		(DnaHmmOrderNumStateMap::Index(1, 3))
		(DnaHmmOrderNumStateMap::Index(2, 3))
		(DnaHmmOrderNumStateMap::Index(3, 3))
		;
	CheckTrainingIncreasesLikelihoods checker(hmm_indices, 3);

	std::for_each(test_seqs.begin(), test_seqs.end(), checker);
}

void
check_hmm_serialization()
{
	cout << "******* check_hmm_serialization()" << endl;

	using namespace boost::archive;

	DnaHmmOrderNumStateMap hmm_map;
	hmm_map.insert_model(1, 0, create_dna_model(1, 0));
	hmm_map.get_model(1, 0);

	//train it on a sequence
	const std::string test_seq = "AGCTATCAGTCGATGATCGATGCTAGTCG";
	hmm_map.get_model(1, 0).train(test_seq);
	BOOST_CHECK_CLOSE(
		hmm_map.get_model(1, 0).get_likelihood(test_seq),
		hmm_map.get_model(1, 0).get_likelihood(test_seq),
		0.001);

	{
		std::stringstream stream;
		text_oarchive(stream) << const_cast<const DnaHmmOrderNumStateMap &>(hmm_map);
		DnaHmmOrderNumStateMap hmm_map_copy;
		text_iarchive(stream) >> hmm_map_copy;

		BOOST_CHECK_CLOSE(
			hmm_map			.get_model(1, 0).get_likelihood(test_seq),
			hmm_map_copy	.get_model(1, 0).get_likelihood(test_seq),
			0.001);
	}

	{
		std::stringstream stream;
		binary_oarchive(stream) << const_cast<const DnaHmmOrderNumStateMap &>(hmm_map);
		DnaHmmOrderNumStateMap hmm_map_copy;
		binary_iarchive(stream) >> hmm_map_copy;

		BOOST_CHECK_CLOSE(
			hmm_map			.get_model(1, 0).get_likelihood(test_seq),
			hmm_map_copy	.get_model(1, 0).get_likelihood(test_seq),
			0.001);
	}
}

void register_hmm_dna_tests(boost::unit_test::test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_hmm_serialization), 0);
	test->add(BOOST_TEST_CASE(&check_hmm_single_training), 0);
	test->add(BOOST_TEST_CASE(&check_unsigned_power), 0);
	test->add(BOOST_TEST_CASE(&check_random_dna_creation), 0);
	test->add(BOOST_TEST_CASE(&check_emission_dna_conversion), 0);
	test->add(BOOST_TEST_CASE(&check_hmm_map), 0);

	//test->add(BOOST_TEST_CASE(&compare_hmm_sizes), 0);
}

