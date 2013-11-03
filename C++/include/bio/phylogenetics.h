

#ifndef BIO_PHYLOGENETICS_H_
#define BIO_PHYLOGENETICS_H_

#include "bio/defs.h"
#include "bio/sequence.h"

BIO_NS_START


/** The details of one phylogenetic sequence. */
struct PhylogeneticSequence
{
	/** How similar is the sequence to the subject sequence? (between 0 and 1). */
	double conservation;

	/** The actual sequence. */
	seq_t sequence;

	PhylogeneticSequence(const seq_t & sequence, double conservation)
		: conservation(conservation)
		, sequence(sequence)
	{
	}
};




BIO_NS_END

#endif //BIO_PHYLOGENETICS_H_

