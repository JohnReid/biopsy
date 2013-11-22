/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
USING_BIO_NS

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/regex.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;





struct ShowReMosApp : Application
{
	std::string remo_extraction_filename;
	std::string sequence;

	ShowReMosApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_space.bin"), "remo extraction file")
			("sequence,s", po::value(&sequence)->default_value("."), "regex to match sequence name")
			;

		get_positional_options().add("sequence", 1);
	}

	int task()
	{
		if ("" == sequence)
		{
			throw std::logic_error( "No regex defined to match gene against" );
		}
		else
		{
			cout << "Showing remos that match regex: " << sequence << "\n";
		}

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

		//show the remos that match the sequence
		unsigned num_matched = 0;
		regex re(sequence);
		smatch what;
		{
			const unsigned num_groups = remo->sequence_groups.size();
			cout << "Searching " << num_groups << " sequence groups\n";

			for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
				remo->sequence_groups.end() != sg;
				++sg)
			{
				bool matched = false;
				//check each of the sequences in the group to see if the id matches
				for (ReMoSequence::list_t::const_iterator s = sg->get()->sequences.begin();
					sg->get()->sequences.end() != s && ! matched;
					++s)
				{
					matched = regex_search(s->get()->id, what, re);
				}

				if (matched)
				{
					//print out group
					cout << **sg << "\n";
					++num_matched;
				}
			}
		}
		cout << "Matched " << num_matched << " groups\n";

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return ShowReMosApp().main(argc, argv);
}

