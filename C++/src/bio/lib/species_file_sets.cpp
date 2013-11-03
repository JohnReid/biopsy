/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/sequence_collection.h"
#include "bio/species_file_sets.h"
#include "bio/environment.h"
#include "bio/random.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <sstream>
#include <ios>
using namespace std;

BIO_NS_START

static species_file_sets_t _species_file_sets;


typedef boost::shared_ptr<boost::filesystem::ifstream> ifstream_ptr;
typedef std::map<boost::filesystem::path, ifstream_ptr> StreamMap;
static StreamMap stream_map;

boost::filesystem::ifstream & get_stream(const boost::filesystem::path & file)
{
	//look for the file in the map...
	StreamMap::iterator stream = stream_map.find(file);
	if (stream == stream_map.end())
	{
		stream = stream_map.insert(StreamMap::value_type(file, ifstream_ptr(new fs::ifstream(file)))).first;
	}
	if (! stream->second->good())
	{
		throw std::string("Could not open file: ") + file._BOOST_FS_NATIVE();
	}
	return *(stream->second);
}


/** Make sure we have built the file sets of chromosomes for each species. */
void
build_species_file_sets()
{
	static bool already_done = false;
	if (! already_done)
	{
		//for each species
		for (std::vector<std::string>::const_iterator i = BioEnvironment::singleton().species_prefixes.begin();
			i != BioEnvironment::singleton().species_prefixes.end();
			++i)
		{
			boost::filesystem::path
				chromo_path(
					BioEnvironment::singleton().get_chromosome_dir()
				);

			_species_file_sets[*i] =
				ChromosomesFileSet::ptr_t(
					new ChromosomesFileSet(
						chromo_path,
						*i));
		}

		already_done = true;
	}
}

species_file_sets_t &
get_species_file_sets()
{
	static bool already_built = false;
	if (! already_built)
	{
		build_species_file_sets();

		already_built = true;
	}

	return _species_file_sets;
}

#ifdef WIN32

std::string
get_last_error()
{
	LPVOID lpMsgBuf;
	DWORD dw = GetLastError(); 

	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		FORMAT_MESSAGE_FROM_SYSTEM,
		NULL,
		dw,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &lpMsgBuf,
		0, NULL );

	stringstream error_msg;
	error_msg
		<< dw << ": "
		<< reinterpret_cast<const char *>(lpMsgBuf);

	LocalFree(lpMsgBuf);

	return error_msg.str();
}

#endif

void
get_random_sequences(unsigned number, unsigned length, SeqList & seq_list)
{
	//for each sequence we want to get
	for (unsigned i = 0; number != i; ++i)
	{
		//choose a species file set at random
		species_file_sets_t::const_iterator file_set = get_species_file_sets().begin();
		for (unsigned species_idx = get_uniform_index(get_species_file_sets().size());
			0 != species_idx;
			--species_idx)
		{
			++file_set;
		}
		//std::cout << "Have chosen " << file_set->first << " file set" << std::endl;
		
		//keep looking until we find a sequence of the given length in it
		seq_t sequence;
		do
		{
			//get a random file from the set
			const boost::filesystem::path file = file_set->second->next_random_file();

			//choose a random position in the file
			file_it random_position;
			{
				//open the file
				file_it file_start(file._BOOST_FS_NATIVE());
				if (! file_start)
				{
#ifdef WIN32

					throw BIO_MAKE_STRING(
						"Could not open file: "
							<< file._BOOST_FS_NATIVE()
							<< ": "
							<< get_last_error());

#else

					throw std::string("Could not open file: ") + file._BOOST_FS_NATIVE();

#endif //WIN32
				}

				//find the end of the file
				file_it file_end = file_start.make_end();

				//skip first line of fasta format
				file_it seq_start = file_start;
				while (*seq_start != '\n') {
					++seq_start;
				}
				++seq_start;

				//get the size of the sequence in the file
				size_t seq_size = size_t(boost::filesystem::file_size(file) - (seq_start - file_start));

				//go to a random position in the sequence
				size_t random_index = get_uniform_index(seq_size);
				random_position = seq_start + random_index;
			}

			//predicate to test whether we have a good base
			is_known_nucleotide is_known;

			//keep adding bases until we get our sequence length or find an unknown
			sequence.clear();
			while (
				(
					is_known(*random_position)
					|| '\r' == *random_position //ignore new lines
					|| '\n' == *random_position
				)
				&& sequence.size() < length)
			{
				if ('\r' != *random_position && '\n' != *random_position)
				{
					sequence.push_back(*random_position);
				}
				++random_position;
			}
		}
		while (sequence.size() < length);
		BOOST_ASSERT(length == sequence.size());

		seq_list.push_back(sequence);
	}

	BOOST_ASSERT(number == seq_list.size());
}

SequenceCollection::ptr_t
get_random_sequence_collection(unsigned number, unsigned desired_length)
{
	std::auto_ptr<SequenceCollectionVector> result(new SequenceCollectionVector());

	//for each sequence we want to get
	for (unsigned i = 0; number != i; ++i)
	{
		//choose a species file set at random
		species_file_sets_t::const_iterator file_set = get_species_file_sets().begin();
		for (unsigned species_idx = get_uniform_index(get_species_file_sets().size());
			0 != species_idx;
			--species_idx)
		{
			++file_set;
		}
		//std::cout << "Have chosen " << file_set->first << " file set" << std::endl;
		
		//keep looking until we find a sequence of the given length in it
		while (true)
		{
			//get a random file from the set
			const boost::filesystem::path file = file_set->second->next_random_file();

			//choose a random position in the file
			unsigned file_size = unsigned(fs::file_size(file));
			unsigned start_position = get_uniform_index(file_size - desired_length);

			//open the file
			fs::ifstream & file_stream = get_stream(file);

			//go to the random position
			file_stream.seekg(start_position, std::ios_base::beg);

			//predicate to test whether we have a good base
			is_known_nucleotide is_known;

			//keep adding bases until we get our sequence length or find an unknown
			seq_t sequence;
			while (desired_length != sequence.size() && file_stream)
			{
				char c = file_stream.get();

				//only allow known characters and new lines
				if (! is_known(c) && '\r' != c && '\n' != c)
				{
					break;
				}

				if (is_known(c))
				{
					sequence.push_back(c);
				}
			}

			//did we find enough known chars?
			if (desired_length == sequence.size())
			{
				result->add_sequence(sequence);
				break;
			}
		}
	}

	BOOST_ASSERT(number == result->num_sequences());

	return SequenceCollection::ptr_t(result.release());
}



BIO_NS_END
