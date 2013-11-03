/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/fasta.h"
#include "FastaLexer.hpp"
#include "NucleoLexer.hpp"
#include "FastaParser.hpp"

#include <boost/algorithm/string/trim.hpp>

using namespace antlr;
using namespace std;


BIO_NS_START


void
parse_fasta(
	istream & stream,
	stringstream & sequence)
{
	FastaLexer fasta_lexer(stream);
	NucleoLexer nucleo_lexer(stream);

	antlr::TokenStreamSelector selector;
	selector.addInputStream(&fasta_lexer, "fasta");
	selector.addInputStream(&nucleo_lexer, "nucleo");
	selector.select("fasta"); //start state

	FastaParser parser(selector);
	parser.sps.selector = &selector;
	parser.sequence = &sequence;

	parser.fasta();

}

std::string
parse_fasta_2(
	std::istream & stream,
	fasta_file_map_t & map)
{
	std::string first_index;

	std::string current_index;
	for (std::string line; getline(stream, line); )
	{
		if ("" == line)
		{
			continue;
		}
		if ('>' == line[0])
		{
			current_index = line.substr(1);
			if ("" == first_index)
			{
				first_index = current_index;
			}

			//do we already have this index in our map?
			if (map.end() != map.find(current_index))
			{
				throw BIO_MAKE_STRING("Already have this index in our map: " << current_index);
			}
		}
		else
		{
			boost::trim(line);
			map[current_index].append(line);
		}
	}

	return first_index;
}


BIO_NS_END
