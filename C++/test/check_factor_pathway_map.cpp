
#include <bio/amigo_pathways.h>
#include <bio/factor.h>
#include <bio/biobase_db.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <string>
#include <iterator>
#include <iostream>
#include <set>
using namespace std;

//#define VERBOSE_CHECKING



void
check_factor_pathway_map()
{
	cout << "******* check_factor_pathway_map(): Checking " << BiobaseDb::singleton().get_factors().size() << " factors" << endl;

	set<DatabaseRef> refs_in_factors;
	set<DatabaseRef> refs_in_pathways;

	typedef map<AmigoPathwayPtr, size_t> AmigoPathwayHitCounts;
	AmigoPathwayHitCounts pathway_hit_counts;

	{
#ifdef VERBOSE_CHECKING
		progress_timer progress;
#endif

		size_t num_factors_matched_to_pathways = 0;
		for (Factor::map_t::const_iterator f_it = BiobaseDb::singleton().get_factors().begin();
			f_it != BiobaseDb::singleton().get_factors().end();
			++f_it)
		{
			bool matched_pathway = false;
			const Factor::ptr_t & factor = f_it->second;
			const DatabaseRefVec & database_refs = factor->database_refs;
			for (DatabaseRefVec::const_iterator db_it = database_refs.begin(); db_it != database_refs.end(); ++db_it)
			{
				refs_in_factors.insert(*db_it);

				AmigoPathwayVec containing_pathways;
				pathways.find_pathways_containing(*db_it, inserter(containing_pathways, containing_pathways.begin()));
				if (containing_pathways.size() > 0) {
					matched_pathway = true;
#ifdef VERBOSE_CHECKING
					cout << factor->get_link() <<  ";" << *db_it << ";";
#endif
					for (AmigoPathwayVec::const_iterator i = containing_pathways.begin(); i != containing_pathways.end(); ++i)
					{

#ifdef VERBOSE_CHECKING
						cout << " " << (*i)->get_name();
#endif

						AmigoPathwayHitCounts::iterator pathway_count = pathway_hit_counts.find(*i);
						if (pathway_hit_counts.end() == pathway_count) {
							pathway_hit_counts[*i] = 1;
						} else {
							pathway_count->second++;
						}
					}

#ifdef VERBOSE_CHECKING
					cout << endl;
#endif

				}
			}
			if (matched_pathway) {
				++num_factors_matched_to_pathways;
			}
		}

#ifdef VERBOSE_CHECKING
		cout << "Found " << refs_in_factors.size() << " unique database refs in factors" << endl;
		cout << num_factors_matched_to_pathways << " factors had matching pathways" << endl;
#endif

	}

	{

#ifdef VERBOSE_CHECKING
		cout << "Building set of database refs in pathways" << endl;
		progress_timer timer;
#endif

		for (AmigoPathwayMap::const_iterator p_it = pathways.pathways.begin(); p_it != pathways.pathways.end(); ++p_it)
		{
			const AmigoPathwayPtr & pathway = p_it->second;
			const DatabaseRefVec & database_refs = pathway->database_refs;
			for (DatabaseRefVec::const_iterator db_it = database_refs.begin(); database_refs.end() != db_it; ++db_it)
			{
				refs_in_pathways.insert(*db_it);
			}
		}

#ifdef VERBOSE_CHECKING
		cout << "Found " << refs_in_pathways.size() << " unique refs in pathways" << endl;
#endif
	}

	set<DatabaseRef> ref_intersection;
	set_intersection(
		refs_in_pathways.begin(),
		refs_in_pathways.end(),
		refs_in_factors.begin(),
		refs_in_factors.end(),
		inserter(ref_intersection, ref_intersection.begin()));

#ifdef VERBOSE_CHECKING
	cout << ref_intersection.size() << " database refs in intersection" << endl;

	for (AmigoPathwayHitCounts::const_iterator i = pathway_hit_counts.begin(); pathway_hit_counts.end() != i; ++i) {
		cout << i->second << " hits in \"" << i->first->get_name() << "\"" << endl;
	}
#endif
}


void register_factor_pathway_map_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_factor_pathway_map), 0);
}


