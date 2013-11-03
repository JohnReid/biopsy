
#include <bio/pssm_motif.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;

#include <iostream>
using namespace std;

//#define VERBOSE_CHECKING
//BOOST_TEST_DONT_PRINT_LOG_VALUE( seq_t::const_iterator );


void
check_pssm_motif_parse()
{
	cout << "******* check_pssm_motif_parse()\n";

	//check some descriptions we expect to parse...
	{
		PssmMotif::ptr_t motif = PssmMotif::parse("PSSM M4 OR PSSM R00004; [1,2] PSSM M44; ");
		BOOST_CHECK_EQUAL(motif->elements.size(), 2u);
	}
	{
		PssmMotif::ptr_t motif = PssmMotif::parse("PSSM M2 OR PSSM R00002;");
		BOOST_CHECK_EQUAL(motif->elements.size(), 1u);
	}
	{
		PssmMotif::ptr_t motif = PssmMotif::parse("PSSM M3 OR PSSM R00003; PSSM M33; ");
		BOOST_CHECK_EQUAL(motif->elements.size(), 1u);
	}
	{
		PssmMotif::ptr_t motif = PssmMotif::parse("PSSM M1; ");
		BOOST_CHECK_EQUAL(motif->elements.size(), 1u);
	}


	//check some descriptions we don't expect to parse...
	{
		const std::vector< std::string > descriptions =
			assign::list_of
				( "" )
				( "GARBAGE" )
				( "PSSM M1; GARBAGE" )
				;

		for (std::vector< std::string >::const_iterator d = descriptions.begin();
			descriptions.end() != d;
			++d)
		{
			bool parsed = true;
			try
			{
				PssmMotif::ptr_t motif = PssmMotif::parse(*d);
			}
			catch (...)
			{
				parsed = false;
			}

			BOOST_CHECK_EQUAL(parsed, false);
		}
	}
}

void
check_pssm_motif_simple_matches()
{
	cout << "******* check_pssm_motif_simple_matches()\n";

	match_result_vec_t matches;
	matches.push_back(
		MatchResults(
			TableLink(MATRIX_DATA, 1), //length 12
			Hit(.99f, 1)));
	matches.push_back(
		MatchResults(
			TableLink(MATRIX_DATA, 2),
			Hit(.99f, 20)));

	{
		PssmMotif::HitVec hits;
		PssmMotif::parse("PSSM M1")->find_in(matches, hits);
		BOOST_CHECK_EQUAL(hits.size(), 1u);
		BOOST_CHECK_EQUAL(hits[0].size(), 1u);
		BOOST_CHECK_EQUAL(hits[0][0].match_result->link, TableLink(MATRIX_DATA, 1));
	}
	{
		PssmMotif::HitVec hits;
		PssmMotif::parse("PSSM M1 OR PSSM M2")->find_in(matches, hits);
		BOOST_CHECK_EQUAL(hits.size(), 2u);
		BOOST_CHECK_EQUAL(hits[0].size(), 1u);
		BOOST_CHECK_EQUAL(hits[1].size(), 1u);
		BOOST_CHECK_EQUAL(hits[0][0].match_result->link, TableLink(MATRIX_DATA, 1));
		BOOST_CHECK_EQUAL(hits[1][0].match_result->link, TableLink(MATRIX_DATA, 2));
	}
	{
		PssmMotif::HitVec hits;
		PssmMotif::parse("PSSM M1; PSSM M2")->find_in(matches, hits);
		BOOST_CHECK_EQUAL(hits.size(), 1u);
		BOOST_CHECK_EQUAL(hits[0].size(), 2u);
		BOOST_CHECK_EQUAL(hits[0][0].match_result->link, TableLink(MATRIX_DATA, 1));
		BOOST_CHECK_EQUAL(hits[0][1].match_result->link, TableLink(MATRIX_DATA, 2));
	}
	{
		PssmMotif::HitVec hits;
		PssmMotif::parse("PSSM M1; [1,20] PSSM M2")->find_in(matches, hits);
		BOOST_CHECK_EQUAL(hits.size(), 1u);
		BOOST_CHECK_EQUAL(hits[0].size(), 2u);
		BOOST_CHECK_EQUAL(hits[0][0].match_result->link, TableLink(MATRIX_DATA, 1));
		BOOST_CHECK_EQUAL(hits[0][1].match_result->link, TableLink(MATRIX_DATA, 2));
	}
	{
		PssmMotif::HitVec hits;
		PssmMotif::parse("PSSM M1; [1,2] PSSM M2")->find_in(matches, hits);
		BOOST_CHECK_EQUAL(hits.size(), 0u);
	}
}

void
register_pssm_motif_tests(boost::unit_test::test_suite * test)
{
	//test->add( BOOST_TEST_CASE( &check_pssm_motif_parse ), 0);
	test->add( BOOST_TEST_CASE( &check_pssm_motif_simple_matches ), 0);
}


