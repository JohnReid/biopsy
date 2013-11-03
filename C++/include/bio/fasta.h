#ifndef BIOBASE_FASTA_H_
#define BIOBASE_FASTA_H_

#include "bio/defs.h"
#include "bio/sequence.h"

#include <sstream>
#include <iostream>


BIO_NS_START

void parse_fasta(
	std::istream & stream,
	std::stringstream & sequence);




/** Maps identifier lines in a fasta file to their sequences. */
typedef std::map< std::string, seq_t > fasta_file_map_t;

/** Suitable for parsing reasonably sized fasta files with one or more sequences in them... Returns key for first sequence in file. */
std::string parse_fasta_2(
	std::istream & stream,
	fasta_file_map_t & map);



BIO_NS_END


#endif //BIOBASE_FASTA_H_
