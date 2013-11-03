/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo_analysis.h"
#include "bio/run_match.h"
#include "bio/biobase_db.h"
#include "bio/biobase_binding_model.h"
#include "bio/biobase_score.h"
#include "bio/adjust_hits.h"
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





struct AnalyseReMoExtractionApp : Application
{
	std::string remo_extraction_filename;
	std::string results_filename;
	float_t threshold;
	float_t phylo_threshold;
	bool masked;
	bool use_consensus_sequences;
	bool output_bifa_analysis;
	string species_filter;
	string matrix_regex;

	AnalyseReMoExtractionApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_space.bin"), "remo extraction file")
			("output,o", po::value(&results_filename)->default_value("remo_analysis.bin"), "results file")
			("threshold,t", po::value(&threshold)->default_value(0.05f), "threshold for hits")
			("phylo_threshold,p", po::value(&phylo_threshold)->default_value(0.05f), "threshold for phylogenetic sequences")
			("species_filter", po::value(&species_filter)->default_value("V"), "species filter")
			("consensus,c", po::value(&use_consensus_sequences)->default_value(true), "use consensus_sequences")
			("matrix_regex,m", po::value(&matrix_regex)->default_value("."), "RegEx to match matrix names")
			("masked", po::value(&masked)->default_value(true), "use masked remos")
			("output_bifa_analysis,b", po::value(&output_bifa_analysis)->default_value(false), "output in bifa analysis format")
			;
	}

	int task()
	{
		cout << "Reading remos from " << remo_extraction_filename << endl;
		cout << "Using threshold of " << threshold << endl;
		cout << "Using phylo threshold of " << phylo_threshold << endl;
		cout << "Analysing " << (masked ? "masked" : "unmasked") << " remos\n";
		cout << "Using species filter: " << species_filter << "\n";
		cout << "Using matrix regex: " << matrix_regex << "\n";
		cout << (use_consensus_sequences ? "Using" : "Not using") << " consensus sequences\n";
		cout << (masked ? "Using" : "Not using") << " masked remos\n";
		cout << "Output will be in " << (output_bifa_analysis ? "bifa" : "remo") << " format\n";


		//deserialise the binary remo archive
		ReMoExtraction::ptr_t remo;
		{
			fs::path remo_extraction_archive( remo_extraction_filename );
			cout << "Deserialising remo extraction from \"" << remo_extraction_archive._BOOST_FS_NATIVE() << "\"\n";
			boost::progress_timer timer;
			remo = ReMoExtraction::deserialise(remo_extraction_archive);
		}


		//create the predicates to choose which matrices and sites we use
		BiobasePssmFilter filter(
			use_consensus_sequences,
			species_filter,
			matrix_regex );

		//do the analysis on each remo
		ReMoAnalysis::map_t analysis_map;
		{
			const unsigned num_groups = remo->sequence_groups.size();
			BiobaseDb::singleton();
			cout << "Analysing remos (" << num_groups << " sequence groups)\n";
			boost::progress_display show_progress( num_groups );
			boost::progress_timer timer;

			for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
				remo->sequence_groups.end() != sg;
				++sg, ++show_progress)
			{
				try
				{
					for (ReMoBundle::map_t::const_iterator rb = sg->get()->remo_bundles.begin();
						sg->get()->remo_bundles.end() != rb;
						++rb)
					{
						//check we don't already have analysis for this remo
						const std::string sequence = rb->first.id;
						const ReMoLocation location = sg->get()->get_sequence_for(sequence)->location;
						const ReMoRange range = rb->first.range;
						if (analysis_map[sequence][location].end() != analysis_map[sequence][location].find(range))
						{
							throw BIO_MAKE_STRING("Already have analysis for: " << sequence << " " << location << " " << range);
						}

						ReMoAnalysis::ptr_t analysis(new ReMoAnalysis);

						//get and score the centre sequence
						analysis->sequence = rb->second->get_sequence(rb->second->centre_sequence, masked);
						score_all_biobase_pssms(
							make_sequence_scorer(
								analysis->sequence.begin(),
								analysis->sequence.end(),
								threshold,
								std::inserter( analysis->results, analysis->results.begin() )
							),
							filter,
							Link2BiobaseBindingModel() 
						);

						//get the set of sequence ids to adjust centre sequence score with
						const ReMoBundle::id_set_t ids = rb->second->get_sequence_ids();

						SeqList sequences;
						for (ReMoBundle::id_set_t::const_iterator id = ids.begin();
							ids.end() != id;
							++id)
						{
							//don't rescore the centre sequence
							if( rb->second->centre_sequence != *id )
							{
								sequences.push_back( rb->second->get_sequence(*id, masked) );
							}
						}

						adjust_hits(
							analysis->results,
							sequences,
							phylo_threshold);

						remove_under_threshold(
							analysis->results,
							phylo_threshold );

						analysis_map[sequence][location][range] = analysis;
					}

				}
				catch (const exception & ex)
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
			cout << "Got analysis for " << analysis_map.size() << " remos\n";
		}

		//serialise
		{
			cout << "Serialising remo analysis to \"" << results_filename << "\"\n";
			std::ofstream stream(results_filename.c_str(), std::ios::binary);
			if (output_bifa_analysis)
			{
				BiFaAnalysis::map_t bifa_map;
				ReMoAnalysis::convert(analysis_map, bifa_map);
				boost::archive::binary_oarchive(stream) << const_cast<const BiFaAnalysis::map_t &>(bifa_map);
			}
			else
			{
				boost::archive::binary_oarchive(stream) << const_cast<const ReMoAnalysis::map_t &>(analysis_map);
			}
		}

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return AnalyseReMoExtractionApp().main(argc, argv);
}

