/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
#include "bio/run_match.h"
#include "bio/biobase_db.h"
#include "bio/counter.h"
#include "bio/markov_model.h"
USING_BIO_NS

#include <boost/progress.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/archive/binary_oarchive.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;





struct ReMoInfoApp : Application
{
	std::string remo_extraction_filename;
	unsigned length_resolution;
	bool masked;
	bool show_slice_info;
	bool show_positions;
	bool show_lengths;
	bool show_conservations;

	ReMoInfoApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_space.bin"), "remo extraction file")
			("masked,m", po::value(&masked)->default_value(true), "use masked remos")
			("show_slice_info,s", po::value(&show_slice_info)->default_value(false), "show slice info")
			("show_positions,p", po::value(&show_positions)->default_value(true), "show positions")
			("show_lengths,l", po::value(&show_lengths)->default_value(true), "show remo lengths")
			("show_conservations,c", po::value(&show_conservations)->default_value(true), "show conservation per species")
			("length_resolution", po::value(&length_resolution)->default_value(100), "resolution of remo length counts")
			;
	}

	int task()
	{
		cout << (masked ? "Using" : "Not using") << " masked remos\n";
		cout << (show_positions ? "Showing" : "Not showing") << " remo positions\n";
		cout << (show_lengths ? "Showing" : "Not showing") << " remo lengths\n";

		//deserialise the binary remo archive
		ReMoExtraction::ptr_t remo;
		{
			fs::path
				remo_extraction_archive(
					remo_extraction_filename
				);
			cout << "Deserialising remo extraction from \"" << remo_extraction_archive._BOOST_FS_NATIVE() << "\"\n";
			boost::progress_timer timer;
			remo = ReMoExtraction::deserialise(remo_extraction_archive);
		}


		typedef Counter< ReMoSpecies > species_counter_t;
		species_counter_t species_counter;

		typedef Counter< unsigned > c_ang_g_content_counter_t;
		c_ang_g_content_counter_t c_ang_g_content_counter;

		typedef Counter< unsigned > unknown_counter_t;
		unknown_counter_t unknown_counter;

		typedef Counter< unsigned > length_counter_t;
		length_counter_t length_counter;

		typedef Counter< unsigned > conservation_counter_t;
		typedef std::map< ReMoSpecies, conservation_counter_t > species_conservation_map_t;
		species_conservation_map_t species_conservations;

		typedef Counter< int > position_counter_t;
		position_counter_t upstream_position_counter;
		position_counter_t downstream_position_counter;
		position_counter_t gene_region_position_counter;

		unsigned num_remos = 0;
		unsigned total_remo_bases = 0;
		unsigned total_num_species = 0;
		const unsigned num_groups = remo->sequence_groups.size();
		std::set< std::string > genes;
		unsigned num_unparsed_genes = 0;

		//do the analysis on each remo
		{

			for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
				remo->sequence_groups.end() != sg;
				++sg)
			{
				//try and parse the centre sequence id
				ReMoSequenceId seq_id;
				if (! parse_sequence_id(sg->get()->get_centre_sequence().id, seq_id))
				{
					++num_unparsed_genes;
				}
				else
				{
					genes.insert(seq_id.gene_id);
				}

				try
				{
					for (ReMoBundle::map_t::const_iterator rb = sg->get()->remo_bundles.begin();
						sg->get()->remo_bundles.end() != rb;
						++rb, ++num_remos)
					{
						//for each species in the remo
						for (ReMo::map_t::const_iterator s = rb->second->remos.begin();
							rb->second->remos.end() != s;
							++s, ++total_num_species)
						{
							const ReMoSpecies species = sg->get()->get_sequence_for(s->first)->species;
							species_conservations[species].increment(s->second.begin()->get()->conservation);
							species_counter.increment(species);
						}

						//get the centre sequence
						const seq_t centre_sequence = rb->second->get_sequence(rb->second->centre_sequence, masked);

						//increment total number of bases
						total_remo_bases += centre_sequence.size();

						//look at the c and g content
						{
							DnaSymbolAlphabet alphabet;
							MarkovModel< 0 > mm(4);
							mm.add_to_counts(centre_sequence.begin(), centre_sequence.end(), alphabet);
							const double c_and_g_content = double(mm.counts[alphabet('c')] + mm.counts[alphabet('g')]) / double(mm.total_count);
							c_ang_g_content_counter.increment(unsigned(c_and_g_content * 20));
						}

						//update the position counts
						const int position = rb->second->get_centre_remo().range.get_centre() / 10000;
						switch (sg->get()->get_centre_sequence().location)
						{
						case REMO_LOC_UPSTREAM: upstream_position_counter.increment(position); break;
						case REMO_LOC_DOWNSTREAM: downstream_position_counter.increment(position); break;
						case REMO_LOC_GENEREGION: gene_region_position_counter.increment(position); break;
						case REMO_LOC_UNDEFINED: break;
						default: throw std::logic_error( "Unknown location type" );
						}

						//update the length counts
						length_counter.increment(centre_sequence.size() / length_resolution);

						//update the unknown counter
						const unsigned num_unknown =
							centre_sequence.size()
							- std::count_if(centre_sequence.begin(), centre_sequence.end(), is_known_nucleotide());
						unknown_counter.increment(20 * num_unknown / centre_sequence.size());
					}

				}
				catch (const std::exception & ex)
				{
					cerr << "Error: " << ex.what() << endl;
				}
				catch (const string & msg)
				{
					cerr << "Error: " << msg << endl;
				}
				catch (const char * msg)
				{
					cerr << "Error: " << msg << endl;
				}
				catch (...)
				{
					cerr << "Undefined error" << endl;
				}

			}
		}

		cout << "# genes, " << genes.size() << "\n";
		cout << "# sequence groups, " << num_groups << "\n";
		cout << "# remos, " << num_remos << "\n";
		cout << "# bases, " << total_remo_bases << "\n";
		cout << "\n";

		cout << "Avg # species per remo, " << double(total_num_species) / double(num_remos) << "\n";
		cout << "\n";

		//the distribution of the positions of the remos
		if (show_positions)
		{
			upstream_position_counter.print(
				false,
				std::cout,
				"Upstream position (10 Kb)",
				true,
				0,
				40,
				true);
			cout << "\n";

			gene_region_position_counter.print(
				false,
				std::cout,
				"Gene region position (10 Kb)",
				true,
				0,
				40,
				true);
			cout << "\n";

			downstream_position_counter.print(
				false,
				std::cout,
				"Downstream position (10 Kb)",
				true,
				0,
				40,
				true);
			cout << "\n";
		}

		//the distribution of the conservations of the remos
		if (show_conservations)
		{
			for (species_conservation_map_t::const_iterator s = species_conservations.begin();
				species_conservations.end() != s;
				++s)
			{
				s->second.print(false, cout, BIO_MAKE_STRING(s->first), false, 0, 60);
				cout << "\n";
			}
		}

		//the distribution of the lengths of the remos
		if (show_lengths)
		{
			cout << "Distribution of remo lengths\n";
			length_counter.print(
				false, 
				std::cout, 
				BIO_MAKE_STRING( "Length (" << length_resolution << " bases)" ),
				true,
				0,
				45,
				true );
			cout << "\n";
		}

		cout << "Distribution of unknown base percentages\n";
		unknown_counter.print(false, std::cout, "Unknown bases (5%)");
		cout << "\n";

		//the number of remos in each species
		cout << "Species across remos\n";
		species_counter.print(true, cout, "Species", true);
		cout << "\n";

		//the distribution of c and g content
		cout << "Distribution of remo C and G content\n";
		c_ang_g_content_counter.print(false, cout, "C and G content (5%)", true);
		cout << "\n";


		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return ReMoInfoApp().main(argc, argv);
}

