/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
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

	ReMoInfoApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_space.bin"), "remo extraction file")
			;
	}

	int task()
	{
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

		typedef std::set< std::string > gene_name_set_t;
		gene_name_set_t gene_names;

		//do the analysis on each remo
		{

			for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
				remo->sequence_groups.end() != sg;
				++sg)
			{
				for (ReMoSequence::list_t::const_iterator s = sg->get()->sequences.begin();
					sg->get()->sequences.end() != s;
					++s)
				{
					if (REMO_SPECIES_MOUSE == s->get()->species)
					{
						gene_names.insert(s->get()->id);
					}
				}
			}
		}

		cout << "\n";
		std::copy(gene_names.begin(), gene_names.end(), std::ostream_iterator< std::string >(std::cout, "\n"));
		cout << "\n";

		cout << "# genes, " << gene_names.size() << "\n";
		cout << "\n";

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return ReMoInfoApp().main(argc, argv);
}

