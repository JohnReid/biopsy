/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
#include "bio/markov_model.h"
USING_BIO_NS

#include <boost/filesystem/path.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace std;





struct MarkovModelsApp : Application
{
	std::string remo_extraction_filename;
	bool masked;
	bool sorted;

	MarkovModelsApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_space.bin"), "remo extraction file")
			("masked,m", po::value(&masked)->default_value(true), "use masked remos")
			("sorted,s", po::value(&sorted)->default_value(true), "sort markov model counts")
			;
	}

	int task()
	{
		cout << "Reading remos from " << remo_extraction_filename << "\n";
		cout << (masked ? "Using" : "Not using") << " masked remos\n";
		cout << (sorted ? "Sorting" : "Not sorting") << " markov model counts\n";

		//deserialise the binary remo archive
		ReMoExtraction::ptr_t remo;
		{
			fs::path
				remo_extraction_archive(
					remo_extraction_filename
				);
			cout << "Deserialising remo extraction from \"" << remo_extraction_archive._BOOST_FS_NATIVE() << "\"\n";
			remo = ReMoExtraction::deserialise(remo_extraction_archive);
		}


		DnaSymbolAlphabet alphabet;
		MarkovModel< 0 > zeroth_order_mm(4);
		MarkovModel< 1 > first_order_mm(4);
		MarkovModel< 2 > second_order_mm(4);
		MarkovModel< 3 > third_order_mm(4);


		//do the analysis on each remo
		{
			const unsigned num_groups = remo->sequence_groups.size();

			for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
				remo->sequence_groups.end() != sg;
				++sg)
			{
				for (ReMoBundle::map_t::const_iterator rb = sg->get()->remo_bundles.begin();
					sg->get()->remo_bundles.end() != rb;
					++rb)
				{
					try
					{
						//get the centre sequence
						const seq_t centre_sequence = rb->second->get_sequence(rb->second->centre_sequence, masked);

						zeroth_order_mm.add_to_counts(
							centre_sequence.begin(),
							centre_sequence.end(),
							alphabet);
						first_order_mm.add_to_counts(
							centre_sequence.begin(),
							centre_sequence.end(),
							alphabet);
						second_order_mm.add_to_counts(
							centre_sequence.begin(),
							centre_sequence.end(),
							alphabet);
						third_order_mm.add_to_counts(
							centre_sequence.begin(),
							centre_sequence.end(),
							alphabet);
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
		}

		std::cout << "\n";
		zeroth_order_mm.print(alphabet, sorted); std::cout << "\n";
		first_order_mm.print(alphabet, sorted); std::cout << "\n";
		//second_order_mm.print(alphabet, sorted); std::cout << "\n";
		//third_order_mm.print(alphabet, sorted); std::cout << "\n";

		return 0;
	}

};

int
main(int argc, char * argv[])
{
	return MarkovModelsApp().main(argc, argv);
}

