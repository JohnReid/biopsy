#ifndef BIO_CHROMOSOMES_H_
#define BIO_CHROMOSOMES_H_

#include <boost/spirit/home/classic/iterator/file_iterator.hpp>
#include "bio/defs.h"

#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>

#include <string>
#include <vector>



BIO_NS_START

/** A collection of FASTA chromosome files. Have random access methods proportional to size. */
struct ChromosomesFileSet
{
	typedef boost::shared_ptr<ChromosomesFileSet> ptr_t;

	std::string name;
	std::vector<boost::filesystem::path> files;
	size_t total_size;

	ChromosomesFileSet(
		const boost::filesystem::path & dir_path,
		const std::string & path_prefix,
		const std::string & suffix= ".fa");

	//returns a file with a likelihood proportional to its size
	const boost::filesystem::path & next_random_file() const;

	size_t next_random_index() const;
};

#if BOOST_VERSION >= 103500
 typedef BOOST_SPIRIT_CLASSIC_NS::file_iterator<> file_it;
#else
 typedef boost::spirit::file_iterator<> file_it;
#endif
typedef std::pair<file_it, file_it> file_it_pair;

file_it_pair
get_random_sequence_in(
	const boost::filesystem::path & file,
	size_t training_seq_size);

/** Gets a good (not too many unknowns) sequence at random from the chromosome file set. */
file_it_pair
get_good_random_sequence_in(
	const ChromosomesFileSet & file_set,
	size_t training_seq_size);

struct FileSeq
{
	FileSeq(file_it_pair it_pair);

	file_it_pair it_pair;

	file_it begin() const;
	file_it end() const;

	boost::reverse_iterator<file_it> rbegin() const;
	boost::reverse_iterator<file_it> rend() const;
};
typedef std::vector<FileSeq> FileSeqVec;

void
get_good_random_sequences_in(
	const ChromosomesFileSet & file_set,
	size_t num_seqs,
	size_t training_seq_size,
	FileSeqVec & seqs);

BIO_NS_END

#endif //BIO_CHROMOSOMES_H_
