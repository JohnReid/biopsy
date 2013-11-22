/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include <bio/application.h>
#include <bio/site_data.h>
#include <bio/remo.h>
USING_BIO_NS;


#include <boost/progress.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
using namespace std;



struct ParseReMosApp : Application
{
	std::string text_filename;
	std::string binary_filename;

	ParseReMosApp()
	{
		get_options().add_options()
			("remo,r", po::value(&text_filename)->default_value("remos.txt"), "text remo extraction file")
			("output,o", po::value(&binary_filename)->default_value("remos.bin"), "binary remo extraction file")
			;
	}

	int task()
	{
		ReMoExtraction::ptr_t remo;
		//parse
		cout << "Parsing remo extraction from \"" << text_filename << "\"" << endl;
		{
			boost::progress_timer timer;
			remo = parse_remo_extraction(fs::path(text_filename));
			cout << "Parsed " << remo->sequence_groups.size() << " sequence groups\n";
		}

		fs::path
			remo_extraction_archive(
				binary_filename.c_str()
			);
		cout << "Writing remo extraction to: \"" << remo_extraction_archive._BOOST_FS_NATIVE() << "\"" << endl;
		{
			boost::progress_timer timer;
			remo->serialise(remo_extraction_archive);
		}

		return 0;
	}
};

int
main(int argc, char * argv [])
{
	return ParseReMosApp().main(argc, argv);
}
