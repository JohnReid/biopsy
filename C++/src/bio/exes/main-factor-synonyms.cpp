/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/equivalent_factors.h"
#include "bio/counter.h"
USING_BIO_NS

using namespace boost;
namespace po = boost::program_options;

using namespace std;




struct FactorSynonymsApp : Application
{
	typedef Counter< unsigned > num_factors_counter_t;

	num_factors_counter_t num_factors_counter;
	EquivalentFactors * factor_partition;
	bool show_synonyms;

	FactorSynonymsApp()
		: factor_partition(0)
	{
		get_options().add_options()
			("show_synonyms,s", po::bool_switch(&show_synonyms)->default_value(false), "show synonyms")
			;
	}

	void add_pssm_factors()
	{
		pssm_set pssms;
		get_all_pssms(pssms);


		//
		// Look at each pssm
		//
		std::set< Factor * > referenced_factors; /**< Those factors that pssms refer to */
		Counter< std::string > species_names; /**< The species these factors are for. */
		Counter< std::string > vertebrate_species_names;
		Counter< std::string > invertebrate_species_names;
		for (pssm_set::const_iterator p = pssms.begin();
			pssms.end() != p;
			++p)
		{
			const FactorLinkList & factors = (*p)->get_factors();

			//for each factor
			for (FactorLinkList::const_iterator f = factors.begin();
				factors.end() != f;
				++f)
			{
				BOOST_FOREACH(std::string & species, f->get()->species)
				{
					species_names.increment(species);
				}

				Factor * factor = BiobaseDb::singleton().get_entry<FACTOR_DATA>(f->get()->link);
				if (0 == factor) //we cannot find all the factors for some reason
				{
					//cout << "No factor for " << f->get()->link << "\n";
					//Factor * factor = BiobaseDb::singleton().get_entry<FACTOR_DATA>(f->get()->link);
					continue;
				}

				referenced_factors.insert(factor);

				factor_partition->add_object(f->get()->link.entry_idx);

				if (factor->is_vertebrate())
				{
					BOOST_FOREACH(std::string & species, f->get()->species)
					{
						vertebrate_species_names.increment(species);
					}
				}
				else
				{
					BOOST_FOREACH(std::string & species, f->get()->species)
					{
						invertebrate_species_names.increment(species);
					}
				}
			}
		}

		//species_names.print(true, cout, "Species", true);
		//vertebrate_species_names.print(true, cout, "Vertebrates", true);
		//invertebrate_species_names.print(true, cout, "Invertebrates", true);

		cout << referenced_factors.size() << " factors are referred to by pssms\n";
	}

	void add_car1_factors()
	{
		factor_partition->add_object(726);
		EquivalentFactors::partition_ptr_t partition =
			factor_partition->find_partition(726);
		if (0 == partition)
		{
			throw std::logic_error( "Cannot find partition for factor 726" );
		}
	}

	void add_pou2f1_factors()
	{
		factor_partition->add_object(641);
		factor_partition->add_object(644);
		factor_partition->add_object(1862);
		factor_partition->add_object(1863);
	}

	void add_all_factors()
	{
		cout << "\nNow adding all the factors\n\n";
		for (Factor::map_t::const_iterator f = BiobaseDb::singleton().get_factors().begin();
			BiobaseDb::singleton().get_factors().end() != f;
			++f)
		{
			factor_partition->add_object(f->second->accession_number.entry_idx);
		}
	}

	void analyse_pssms()
	{
		pssm_set pssms;
		get_all_pssms(pssms);

		//
		// Look at each pssm
		//
		unsigned unfound_factors = 0;
		std::set< Factor * > referenced_factors; /**< Those factors that pssms refer to */
		for (pssm_set::const_iterator p = pssms.begin();
			pssms.end() != p;
			++p)
		{
			std::set< EquivalentFactors::partition_ptr_t > factors_for_this_pssm;

			const FactorLinkList & factors = (*p)->get_factors();

			//for each factor
			for (FactorLinkList::const_iterator f = factors.begin();
				factors.end() != f;
				++f)
			{
				EquivalentFactors::partition_ptr_t partition = factor_partition->find_partition((*f)->link.entry_idx);
				if (0 == partition)
				{
					++unfound_factors;
					continue;
					throw BIO_MAKE_STRING("Could not find partition for " << (*f)->link);
				}

				referenced_factors.insert(BiobaseDb::singleton().get_entry<FACTOR_DATA>((*f)->link));

				factors_for_this_pssm.insert(partition);
			}

			num_factors_counter.increment(factors_for_this_pssm.size());

			if (factors_for_this_pssm.size() > 10)
			{
				cout << (*p)->get_link() << " has " << factors_for_this_pssm.size() << " distinct factors\n";
			}
		}

		cout << "Could not find partitions for " << unfound_factors << " factor links\n";

		cout << pssms.size() << " pssms reference " << referenced_factors.size() << " factors\n";
	}


	int task()
	{
		BiobaseDb::singleton();

		factor_partition = &EquivalentFactors::singleton();

		cout << "Have " << BiobaseDb::singleton().get_factors().size() << " factors in Biobase\n";

		analyse_pssms();

#if TESTING
		add_car1_factors();
		print_stats();

		add_pou2f1_factors();
		print_stats();
		sanity_check();

		add_pssm_factors();
		print_stats();
		sanity_check();

		add_all_factors();
		print_stats();
		sanity_check();

#endif

		if (show_synonyms)
		{
			print_partitions();
		}

		print_stats();

		return 0;
	}

	void print_partitions()
	{
		//for each partition
		for (EquivalentFactors::partition_set_iterator p = factor_partition->begin();
			factor_partition->end() != p;
			++p)
		{
			set< string > names;
			for (EquivalentFactors::iterator f = p->get()->begin();
				p->get()->end() != f;
				++f)
			{
				//insert the factor's name and synonyms
				Factor * factor = FactorEquivalence::get_factor(*f);
				names.insert(factor->get_name());
				copy(factor->synonyms.begin(), factor->synonyms.end(), inserter(names, names.begin()));
			}

			//print the names
			copy(names.begin(), names.end(), ostream_iterator< string >(cout, ", "));
			cout << "\n\n";
		}
	}

	void print_stats()
	{
		cout << "Found " << factor_partition->num_partitions() << " factor equivalence class(es)\n";
		EquivalentFactors::partition_ptr_t pou2f1 = factor_partition->find_partition(641);
		if (0 != pou2f1)
		{
			cout << "pou2f1 has " << pou2f1->size() << " equivalent factors\n";
			for (EquivalentFactors::partition_t::const_iterator f = pou2f1->begin();
				pou2f1->end() != f;
				++f)
			{
				Factor * factor = FactorEquivalence::get_factor(*f);
				cout << factor->get_name() << "\n";
			}
		}
		cout << "\n";

		num_factors_counter.print(
			false,
			std::cout,
			"# factors",
			true,
			0,
			50);
	}

	void sanity_check()
	{
		//check pou2f has only one equivalence class
		if (1 != factor_partition->which_partitions(641).size())
		{
			throw std::logic_error( "Bad number of equivalence classes for pou2f1" );
		}
	}
};

int
main(int argc, char * argv[])
{
	return FactorSynonymsApp().main(argc, argv);
}

