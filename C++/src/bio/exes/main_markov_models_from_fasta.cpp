/**
@file

Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/markov_model.h"
#include "bio/fasta.h"
USING_BIO_NS

#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
using namespace boost;
using namespace BOOST_SPIRIT_CLASSIC_NS;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <vector>
#include <strstream>
#include <fstream>
using namespace std;


//typedef BOOST_SPIRIT_CLASSIC_NS::file_iterator<> file_it;
typedef file_iterator<> file_it;


/**
 * Learns Markov models from the sequences in a FASTA file.
 */
struct FastaMarkovModelsApp : Application
{
	typedef std::vector< std::string > filename_vec_t;

	bool sorted;
	filename_vec_t input_filenames;

	FastaMarkovModelsApp()
	{
		get_options().add_options()
			("sorted", po::value(&sorted)->default_value(true), "sort markov model counts")
		    ("input-file", po::value(&input_filenames), "input file")
			;

		get_positional_options().add("input-file", -1);
	}

	int task()
	{

		DnaSymbolAlphabet alphabet;
		MarkovModel< 0 > zeroth_order_mm(4);
		MarkovModel< 1 > first_order_mm(4);
		MarkovModel< 2 > second_order_mm(4);
		MarkovModel< 3 > third_order_mm(4);

		for (filename_vec_t::const_iterator f = input_filenames.begin();
			input_filenames.end() != f;
			++f)
		{
			//open the file
			file_it file_start(*f);
			if (! file_start) {
				throw BIO_MAKE_STRING("Could not open " << *f);
			}

			//find the end of the file
			file_it file_end = file_start.make_end();

			//skip first line of fasta format
			file_it seq_start = file_start;
			while (*seq_start != '\n') {
				++seq_start;
			}
			++seq_start;

			std::cout << *f << "\n";

			zeroth_order_mm.add_to_counts(
				file_start,
				file_end,
				alphabet);
			first_order_mm.add_to_counts(
				file_start,
				file_end,
				alphabet);
#if 0
			second_order_mm.add_to_counts(
				file_start,
				file_end,
				alphabet);
			third_order_mm.add_to_counts(
				file_start,
				file_end,
				alphabet);
#endif

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
	return FastaMarkovModelsApp().main(argc, argv);
}

