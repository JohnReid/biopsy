/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/tss_data.h"

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/tokenizer.hpp>
namespace fs = boost::filesystem;
using namespace boost;

#include <iostream>
using namespace std;


BIO_NS_START

namespace
{
	int parse_int(const string & text)
	{
		stringstream sstream;
		sstream << text;
		int result;
		sstream >> result;
		if (! sstream)
		{
			throw BIO_MAKE_STRING("Could not parse: \"" << text << "\"");
		}
		return result;
	}
};

void TSS::parse_files(
	boost::filesystem::path tss_file,
	boost::filesystem::path clone_file,
	TSS::map_t & map)
{
	map.clear();

	typedef char_separator< char > separator;
	typedef tokenizer< separator > tokenizer;
	separator sep(";");

	//read the clone file...
	{
		const fs::path & file = clone_file;
		if (! fs::exists(file))
		{
			throw BIO_MAKE_STRING(file._BOOST_FS_NATIVE() << " file does not exist");
		}
		fs::ifstream stream(file);
		unsigned line_num = 0;
		while (stream && ! stream.eof())
		{
			string line;
			getline(stream, line);
			if (! stream || stream.eof())
			{
				break;
			}
			++line_num;
			tokenizer tokens(line, sep);
			tokenizer::iterator t = tokens.begin();

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string gene = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string transcript = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string alt_start_num = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string source = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string clone_id = *t++;

			//get a reference to the TSS
			TSS & tss = map[gene];

			//resize the clone vector if necessary
			while (tss.clones.size() < unsigned(parse_int(alt_start_num)))
			{
				tss.clones.push_back(Clone());
			}

			//get a reference to the clone
			Clone & clone = tss.clones[parse_int(alt_start_num) - 1];
			clone.transcript = transcript;

			//update the clone
			if ("undef" == source)
			{
				clone.type = NO_CLONE;
			}
			else if ("RIKEN" == source)
			{
				clone.type = RIKEN_CLONE;
				clone.id = clone_id;
			}
			else if ("TIGR" == source)
			{
				clone.type = TIGR_CLONE;
				clone.id = clone_id;
			}
			else
			{
				throw BIO_MAKE_STRING(source << ": bad clone type with id: " << clone_id);
			}
		}
	}


	//read the tss file...
	{
		const fs::path & file = tss_file;
		if (! fs::exists(file))
		{
			throw BIO_MAKE_STRING(file._BOOST_FS_NATIVE() << " file does not exist");
		}
		fs::ifstream stream(file);
		unsigned line_num = 0;
		while (stream && ! stream.eof())
		{
			string line;
			getline(stream, line);
			if (! stream || stream.eof())
			{
				break;
			}
			++line_num;
			tokenizer tokens(line, sep);
			tokenizer::iterator t = tokens.begin();

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string gene = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string transcript = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string chromosome = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string strand = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string rel_tss_pos = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string tss_pos = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string alt_start_num = *t++;

			if (tokens.end() == t)
				throw BIO_MAKE_STRING("Bad format on line " << line_num << " of " << file._BOOST_FS_NATIVE());
			const string num_start_pos = *t++;

			//get a reference to the TSS
			TSS & tss = map[gene];

			BOOST_ASSERT(tss.clones.size() == unsigned(parse_int(num_start_pos)));

			//get a reference to the clone
			Clone & clone = tss.clones[parse_int(alt_start_num) - 1];
			BOOST_ASSERT(clone.transcript == transcript);
			if ("no_clones" == rel_tss_pos)
			{
				BOOST_ASSERT(NO_CLONE == clone.type);
				BOOST_ASSERT("undef" == tss_pos);
			}
			else if ("no_conclusion" == rel_tss_pos)
			{
				BOOST_ASSERT(NO_CLONE == clone.type);
				clone.type = NO_CONCLUSION_CLONE;
				BOOST_ASSERT("undef" == tss_pos);
			}
			else
			{
				clone.chromosome = chromosome;
				clone.strand = ("1" == strand);
				clone.rel_tss_pos = parse_int(rel_tss_pos);
				clone.tss_pos = parse_int(tss_pos);
			}
		}
	}


	fs::ifstream clone_stream(clone_file);
}

BIO_NS_END



