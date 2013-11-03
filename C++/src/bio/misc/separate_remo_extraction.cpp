/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo_analysis.h"
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





struct SeparateReMoExtractionApp : Application
{
	std::string remo_extraction_filename;
	std::string results_prefix;
	unsigned num_parts;

	SeparateReMoExtractionApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_extraction.bin"), "remo extraction file")
			("num_parts,n", po::value(&num_parts)->default_value(2), "number of parts to split extraction into")
			("results_prefix,p", po::value(&results_prefix)->default_value("remo_extraction_part"), "prefix for results files")
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

		//the vector that holds each part
		typedef std::vector<ReMoExtraction::ptr_t> extraction_vec_t;
		extraction_vec_t extraction_parts;
		for (unsigned i = 0; num_parts != i; ++i)
		{
			extraction_parts.push_back(ReMoExtraction::ptr_t(new ReMoExtraction));
		}

		//iterate through the extraction
		{
			cout << "Splitting extraction\n";
			boost::progress_display show_progress( remo->sequence_groups.size() );
			boost::progress_timer timer;

			unsigned part_idx = 0;
			for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
				remo->sequence_groups.end() != sg;
				++sg, ++show_progress, ++part_idx)
			{
				if (num_parts == part_idx)
				{
					part_idx = 0;
				}

				extraction_parts[part_idx]->sequence_groups.push_back(*sg);
			}
		}

		//serialise
		for (unsigned i = 0; num_parts != i; ++i)
		{
			const std::string output_filename = BIO_MAKE_STRING(results_prefix << '_' << i << ".bin");
			cout << "Serialising remo extraction part to \"" << output_filename << "\"\n";
			std::ofstream stream(output_filename.c_str(), std::ios::binary);
			boost::archive::binary_oarchive(stream) << const_cast<const ReMoExtraction &>(*(extraction_parts[i]));
		}

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return SeparateReMoExtractionApp().main(argc, argv);
}

