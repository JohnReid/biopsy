#include <bio/remo.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <string>
#include <iostream>
using namespace std;

//#define VERBOSE_CHECKING


BOOST_TEST_DONT_PRINT_LOG_VALUE(ReMoExtraction);


void check_remo_parse()
{
	cout << "******* check_remo_parse()" << endl;

	ReMoExtraction::ptr_t remo_extraction;

	//parsing
	{
#ifdef VERBOSE_CHECKING
		cout << "Parsing remos from C:\\data\\ReMos\\CompleteExtraction50.txt\n";
		progress_timer timer;
#endif

		remo_extraction =
			parse_remo_extraction(
				boost::filesystem::path(
					"C:\\data\\ReMos\\CompleteExtraction50.txt"
				)
		);
	}

	//serialization
	{
#ifdef VERBOSE_CHECKING
		cout << "Serialising remos\n";
		progress_timer timer;
#endif

		boost::filesystem::path
			remo_extraction_archive(
				"remo_extraction.bin"
			);

		remo_extraction->serialise(remo_extraction_archive);
		ReMoExtraction::ptr_t remo_extraction_copy = ReMoExtraction::deserialise(remo_extraction_archive);
	}

	//map the genes to remos...
	{
#ifdef VERBOSE_CHECKING
		cout << "Mapping genes to remos\n";
		progress_timer timer;
#endif

		ReMoExtraction::GeneReMoMap gene_remo_map;
		remo_extraction->build_gene_remo_map(gene_remo_map);
	}

#ifdef VERBOSE_CHECKING

	for (ReMoSequenceGroup::list_t::const_iterator sg = remo_extraction->sequence_groups.begin();
		remo_extraction->sequence_groups.end() != sg;
		++sg)
	{
		for (ReMoBundle::map_t::const_iterator rb = sg->get()->remo_bundles.begin();
			sg->get()->remo_bundles.end() != rb;
			++rb)
		{
			ReMoBundle::id_set_t ids = rb->second->get_sequence_ids();
			for (ReMoBundle::id_set_t::const_iterator id = ids.begin();
				ids.end() != id;
				++id)
			{
				cout
					<< *id << ": "
					<< (rb->second->get_centre_sequence_id().id == *id ? "C: " : "P: ")
					<< rb->second->get_sequence(*id, false) << "\n";
			}
			cout << "\n";
		}
		cout << "\n";
	}

#endif
}

void
register_remo_parse_tests(boost::unit_test::test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_remo_parse), 0);
}
