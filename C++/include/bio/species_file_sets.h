
#include "bio/defs.h"
#include "bio/chromosomes_file_set.h"
#include "bio/sequence.h"
#include "bio/sequence_collection.h"

#ifndef BIO_SPECIES_FILE_SETS_H_
#define BIO_SPECIES_FILE_SETS_H_



BIO_NS_START


/** A map from species names to file sets. */
typedef std::map<std::string, ChromosomesFileSet::ptr_t> species_file_sets_t;

/** A map from species names to file sets. */
species_file_sets_t & get_species_file_sets();

/** Get some sequences randomly chosen from the chromosomes. */
void get_random_sequences(unsigned number, unsigned length, SeqList & seq_list);

SequenceCollection::ptr_t
get_random_sequence_collection(unsigned number, unsigned desired_length);

BIO_NS_END


#endif //BIO_SPECIES_FILE_SETS_H_

