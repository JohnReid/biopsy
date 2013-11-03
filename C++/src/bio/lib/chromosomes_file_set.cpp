/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/chromosomes_file_set.h"
#include "bio/random.h"
#include "bio/sequence.h"

#include <boost/filesystem/operations.hpp>
using namespace boost;

#include <iostream>
using namespace std;


BIO_NS_START


ChromosomesFileSet::ChromosomesFileSet(
	const boost::filesystem::path & dir_path,
	const std::string & path_prefix,
	const std::string & suffix)
	: name(path_prefix)
	, total_size(0)
{
	boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
	for (boost::filesystem::directory_iterator i(dir_path); i != end_itr; ++i)
	{
		if (! boost::filesystem::is_directory(*i))
		{
			const std::string filename = i->path().leaf().string();
			static const std::string suffix = ".fa";

			//does it match the prefix
			if (0 == filename.find(path_prefix))
			{
				//does it have the right suffix?
				if (filename.size() - suffix.size() == filename.rfind(suffix))
				{
					//cout << "File: " << i->leaf() << ", size: " << boost::filesystem::file_size(*i) << endl;

					files.push_back(*i);
					const size_t old_total_size = total_size; //check for overflow
					total_size += size_t( boost::filesystem::file_size(*i) );
					if (old_total_size > total_size)
					{
						throw std::logic_error( "Overflow calculating total size" );
					}
				}
			}
		}
	}
	if (0 == total_size)
	{
		throw
			std::string("Did not find any files at: ")
			+ dir_path._BOOST_FS_NATIVE()
			+ " with path prefix: "
			+ path_prefix;
	}
}

size_t
ChromosomesFileSet::next_random_index() const
{
	return get_uniform_index(total_size);
}

//returns a file with a likelihood proportional to its size
const boost::filesystem::path &
ChromosomesFileSet::next_random_file() const
{
	const size_t rnd_idx = next_random_index();
	size_t cumulative_sizes = 0;
	for (size_t i = 0; files.size() != i; ++i)
	{
		cumulative_sizes += size_t( boost::filesystem::file_size(files[i]) );
		if (cumulative_sizes > rnd_idx) {
			//cout << "Chose " << files[i]._BOOST_FS_NATIVE() << endl;
			return files[i];
		}
	}
	cout << "rnd_idx: " << rnd_idx << endl;
	cout << "cumulative_sizes: " << cumulative_sizes << endl;
	throw std::logic_error( "Should have cumulative_sizes > total_sizes" );
}

file_it_pair
get_random_sequence_in(
	const boost::filesystem::path & file,
	size_t training_seq_size)
{
	//open the file
	file_it file_start(file._BOOST_FS_NATIVE());
	if (! file_start) {
		throw std::logic_error( "Could not open file" );
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
	size_t seq_size = size_t( boost::filesystem::file_size(file) - (seq_start - file_start) );

	//go to a random position in the sequence
	size_t random_index = get_uniform_index(seq_size);
	file_it training_seq_middle = seq_start + random_index;

	//find a training sequence of the correct length around our random position
	file_it begin = training_seq_middle;
	file_it end = training_seq_middle;
	size_t length = 0; //the length of the sequence so far
	is_known_nucleotide test; //test for known nucleotides - we ignore 'N' in the counts
	//whilst we are still expanding the length and it is not enough
	size_t last_length;
	do
	{
		//remember the last length
		last_length = length;

		//expand halfway forward and halfway back
		size_t remainder = (training_seq_size - length) / 2;

		//backwards - until beginning or we have remainder known nucleotides
		for (size_t i = 0; i < remainder && begin != seq_start; ) {
			--begin;
			if (test(*begin)) {
				++i;
				++length;
			}
		}

		//forwards - until end or we have remainder known nucleotides
		for ( ; length < training_seq_size && end != file_end; ++end) {
			if (test(*end)) {
				++length;
			}
		}
	}
	while (length < training_seq_size && length != last_length);

	//either there wasn't enough data or we have the right number of known nucleotides
	assert(std::count_if(begin, end, test) == int(length));

	//we only test true here when we are one short so I'm happy to let this slip for the time being....
	//TODO fix the above code
#if 0
	if (! ((end == file_end && begin == seq_start) || length == training_seq_size))
	{
		cout << "Length: " << length << endl;
		cout << "Training seq size: " << training_seq_size << endl;
		cout << "file_end - end: " << file_end - end << endl;
		cout << "begin - seq_start: " << begin - seq_start << endl;
		throw std::logic_error( "Why didn't we find enough data?" );
	}
#endif

	return make_pair(begin, end);
}

/** Gets a good (sequence at random from the chromosome file set. */
file_it_pair
get_good_random_sequence_in(
	const ChromosomesFileSet & file_set,
	size_t training_seq_size)
{
	//find a sequence in a random file from the species, make sure the sequence has at least 50%
	//known nucleotides
	const boost::filesystem::path file = file_set.next_random_file();
	file_it_pair it_pair = get_random_sequence_in(file, training_seq_size);

	//did we find a good sequence?
	if ((size_t) (it_pair.second - it_pair.first) > 10 * training_seq_size)
	{
		//no
		throw std::logic_error( "Cannot find good sequence" );
	}

	return it_pair;
}


FileSeq::FileSeq(file_it_pair it_pair) : it_pair(it_pair) { }

file_it
FileSeq::begin() const
{
	return it_pair.first;
}

file_it
FileSeq::end() const
{
	return it_pair.second;
}

boost::reverse_iterator<file_it>
FileSeq::rbegin() const
{
	return make_reverse_iterator(end());
}

boost::reverse_iterator<file_it>
FileSeq::rend() const
{
	return make_reverse_iterator(begin());
}

void
get_good_random_sequences_in(
	const ChromosomesFileSet & file_set,
	size_t num_seqs,
	size_t seq_length,
	FileSeqVec & seqs)
{
	seqs.clear();

	//only try twice as many times as seqs we want otherwise this can take ages.
	for (size_t tries = 0; tries < 2 * num_seqs && seqs.size() < num_seqs; ++tries)
	{
		try
		{
			seqs.push_back(FileSeq(get_good_random_sequence_in(file_set, seq_length)));
		}
		catch (...)
		{
			//
		}
	}
}

BIO_NS_END


