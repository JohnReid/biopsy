/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#include "bio_test_defs.h"
#include "bio_test_data.h"

#include <bio/biobase_likelihoods.h>
#include <bio/pssm_match.h>
#include <bio/biobase_db.h>
#include <bio/biobase_filter.h>
#include <bio/biobase_data_traits.h>
#include <bio/serialisable.h>
USING_BIO_NS;


#include <boost/progress.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <fstream>
#include <vector>
using namespace std;

#ifdef _DEBUG
#define TEST_SEQ_LENGTH 100
#else
#define TEST_SEQ_LENGTH 10000
#endif


//#define VERBOSE_CHECKING

BOOST_TEST_DONT_PRINT_LOG_VALUE( LikelihoodsCache )


void check_likelihoods_cache()
{
	cout << "******* check_likelihoods_cache()" << endl;

	//create a random sequence
	const size_t test_seq_length = TEST_SEQ_LENGTH;
	seq_t test_seq;
	test_seq.reserve(test_seq_length);
	generate_random_nucleotide_seq(inserter(test_seq, test_seq.begin()), test_seq_length);

	LikelihoodsCache cache;
	{
#ifdef VERBOSE_CHECKING
		cout << "Updating counts" << endl;
		progress_timer timer;
#endif //VERBOSE_CHECKING

		BiobaseDb::singleton();

		cache.update_counts(BiobaseDb::singleton(), test_seq);
	}

	{
#ifdef VERBOSE_CHECKING
		cout << "Serialising cache" << endl;
		progress_timer timer;
#endif //VERBOSE_CHECKING

		serialise< false >( cache, "cache.txt" );
	}

	LikelihoodsCache copy_of_cache;
	{
#ifdef VERBOSE_CHECKING
		cout << "Deserialising cache" << endl;
		progress_timer timer;
#endif //VERBOSE_CHECKING

		deserialise< false >( copy_of_cache, "cache.txt" );
	}
	BOOST_CHECK_EQUAL(cache, copy_of_cache);

	boost::io::ios_precision_saver ips(cout);
	cout.precision(3);

	boost::io::ios_width_saver iws(cout);
	cout.width(8);

	const vector< TableLink > pssms = boost::assign::list_of
		( TableLink( SITE_DATA, 2195 ) )
		( TableLink( MATRIX_DATA, 201 ) )
		( TableLink( MATRIX_DATA, 128 ) )
		( TableLink( MATRIX_DATA, 209 ) )
		//("AG$CONS_02")
		;

	for (vector< TableLink >::const_iterator n = pssms.begin(); pssms.end() != n; ++n)
	{
		const BiobaseCounts * counts = cache.get_counts(*n);

#ifdef VERBOSE_CHECKING
		copy(counts->begin(), counts->end(), ostream_iterator<size_t>(cout, ", ")); cout << endl;
#endif //VERBOSE_CHECKING

		for (int i = 0; i < 2; ++i)
		{
			const bool background = i > 0;
			for (int j = 0; j < 2; ++j)
			{
				const bool or_better = j > 0;
				const BiobaseLikelihoods * likelihoods = cache.get_score_likelihoods(*n, background, or_better);

#ifdef VERBOSE_CHECKING
				cout
					<< *n << " "
					<< (background ? "background " : "binding ")
					<< (or_better ? "or better " : "exact ")
					<< "likelihoods"
					<< endl;
				copy(likelihoods->begin(), likelihoods->end(), ostream_iterator<BIO_NS::float_t>(cout, ", ")); cout << endl;
				cout << endl;
#endif //VERBOSE_CHECKING
			}
		}
	}
}


void check_or_better_likelihoods_bug()
{
	cout << "******* check_or_better_likelihoods_bug()" << endl;

	//const seq_t seq =
	//	"GCTCTCGCTCTTTTTTTTTTTCGCAAAAGGAGGGGAGAGGGGGTAAAAAAATGCTGCACTGTGCGGCGAAGCCGGTGAGTGAGCGGCGC"
	//	"GGGGCCAATCAGCGTGCGCCGTTCCGAAAGTTGCCTTTTATGGCTCGAGCGGCCGCGGCGGCGCCCTATAAAACCCAGCGGCGCGACGCGC"
	//	;

	const seq_t seq =
		"CCCTATAAAACCCAGCG"
		;

	Scorers scorers( seq.begin(), seq.end() );
	hit_vec_t hits;
	Matrix * pssm_M00320 = BiobaseDb::singleton().get_entry< MATRIX_DATA >( 320 );

	scorers.bayesian_or_better_calc(
		*pssm_M00320, 
		std::back_inserter( hits ) );
}

struct check_likelihoods
{
	void operator()( BIO_NS::float_t likelihood ) const
	{
		BOOST_CHECK( likelihood >= BIO_NS::float_t( 0.0 ) );
		BOOST_CHECK( likelihood <= BIO_NS::float_t( 1.0 ) );

		if ( likelihood > BIO_NS::float_t( 1.0 ) )
		{
			cout << "Problem\n";
		}
	}

	void operator()( const BiobaseLikelihoods * likelihoods ) const
	{
		if ( 0 != likelihoods )
		{
			std::for_each(
				likelihoods->begin(),
				likelihoods->end(),
				*this );
		}
	}

	void operator()( const TableLink & pssm ) const
	{
		( *this )( 
			LikelihoodsCache::singleton().get_score_likelihoods(
				pssm,
				false,
				false) );

		( *this )( 
			LikelihoodsCache::singleton().get_score_likelihoods(
				pssm,
				false,
				true) );

		( *this )( 
			LikelihoodsCache::singleton().get_score_likelihoods(
				pssm,
				true,
				false) );

		( *this )( 
			LikelihoodsCache::singleton().get_score_likelihoods(
				pssm,
				true,
				true) );
	}

	template< typename MapValue >
	void operator()( const MapValue & pssm ) const
	{
		( *this )( pssm.first );
	}
};

void check_all_likelihoods()
{
	cout << "******* check_all_likelihoods()" << endl;

	//check_likelihoods()( std::string( "HOX13$CONS_02" ) ); //problem with pssm_likelihoods

	//check each site
	std::for_each(
		get_sites_begin(),
		get_sites_end(),
		check_likelihoods() );

	//check each matrix
	std::for_each(
		get_matrices_begin(),
		get_matrices_end(),
		check_likelihoods() );

}


void register_biobase_likelihoods_tests(test_suite * test)
{
	ensure_biobase_params_built();

	test->add(BOOST_TEST_CASE(&check_all_likelihoods), 0);
	test->add(BOOST_TEST_CASE(&check_likelihoods_cache), 0);
	test->add(BOOST_TEST_CASE(&check_or_better_likelihoods_bug), 0);
}


//
// Are we going to compile this test into its own executable?
//
#ifdef BIO_STANDALONE_TEST
BIO_DEFINE_STANDALONE_TEST( "biobase likeilhoods", register_biobase_likelihoods_tests )
#endif //BIO_STANDALONE_TEST
