
#ifndef BIO_ADJUST_HITS_H_
#define BIO_ADJUST_HITS_H_

#include "bio/defs.h"
#include "bio/binding_model.h"

#include <vector>

BIO_NS_START


/**
Adjust the hits for one phylogenetic sequence.
*/
void
adjust_hits(
	BindingModel::hit_set_t & hits,
	seq_t::const_iterator phylo_begin,
	seq_t::const_iterator phylo_end,
	double threshold,
	BindingModelContext * context = 0 );



typedef std::pair< seq_t::const_iterator, seq_t::const_iterator > seq_it_pair_t;
typedef std::vector< seq_it_pair_t > seq_vec_t;


/**
Adjust the hits for several phylogenetic sequences.
*/
void
adjust_hits(
	BindingModel::hit_set_t & hits,
	const SeqList & sequences,
	double threshold,
	BindingModelContext * context = 0 );


/**
Adjust the hits for several phylogenetic sequences.
*/
void
adjust_hits(
	BindingModel::hit_set_t & hits,
	const std::vector< seq_t > & sequences,
	double threshold,
	BindingModelContext * context = 0 );


/**
Removes all hits under the threshold.
*/
void
remove_under_threshold(
	BindingModel::hit_set_t & hits,
	double threshold );



BIO_NS_END

#endif //BIO_ADJUST_HITS_H_
