/**
@file

Copyright John Reid 2006-2010

*/

#ifndef BIOPSY_ANALYSE_H_
#define BIOPSY_ANALYSE_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"
#include "biopsy/binding_hits.h"

#include "bio/singleton.h"
#include "bio/sequence.h"



#ifndef BIOPSY_ANALYSE_THRESHOLD_DEFAULT
# define BIOPSY_ANALYSE_THRESHOLD_DEFAULT 0.03
#endif //BIOPSY_ANALYSE_THRESHOLD_DEFAULT

#ifndef BIOPSY_BIOBASE_SCORE_THRESHOLD_DEFAULT
# define BIOPSY_BIOBASE_SCORE_THRESHOLD_DEFAULT 0.01
#endif //BIOPSY_ANALYSE_THRESHOLD_DEFAULT




namespace biopsy
{

/**
Scores the pssm on the sequence and returns estimate that the pssm binds in at least one position.
*/
double
score_pssm_on_sequence(
	const std::string & pssm_name,
	const sequence & seq,
	double threshold,
	binding_hit::vec_ptr result );

/**
Generates the biobase scores for the pssm on the sequence.
*/
void
biobase_score_pssm_on_sequence(
	const std::string & pssm_name,
	const sequence & seq,
	double threshold,
	binding_hit::vec_ptr result );

/**
Score a sequence.
*/
binding_hit::vec_ptr
score_pssms_on_sequence(
	const string_vec_ptr & pssm_names,
	const sequence & seq,
	double threshold = BIOPSY_ANALYSE_THRESHOLD_DEFAULT );


/**
Score a sequence returning the biobase scores
*/
binding_hit::vec_ptr
biobase_score_pssms_on_sequence(
	const string_vec_ptr & pssm_names,
	const sequence & seq,
	double threshold = BIOPSY_BIOBASE_SCORE_THRESHOLD_DEFAULT );



/**
0: The adjusted hits for the first sequence.
1: The maximal chain across all sequences.
2: The unadjusted hits for each sequence.
*/
typedef boost::tuple< binding_hit::vec_ptr, binding_hit::vec_ptr, binding_hits_vec_ptr > phylo_sequences_result;


/**
Score the pssms on the sequences.
*/
phylo_sequences_result
score_pssms_on_phylo_sequences(
	string_vec_ptr pssm_names,
	sequence_vec_ptr sequences,
	double threshold,
	double phylo_threshold,
	bool calculate_maximal_chain = true
);


/**
Get the maximum # of sequences for which the phylogenetic analysis algorithm
will bother finding the maximal chain.
*/
unsigned
get_max_chain_max_num_sequences();

/**
Analyse a sequence.
*/
binding_hit::vec_ptr
analyse(
	const sequence & seq,
	double threshold = BIOPSY_ANALYSE_THRESHOLD_DEFAULT );


/**
Analyse a sequence and adjust using a phlyogenetic comparison.
*/
binding_hit::vec_ptr
analyse_phylo(
	const sequence & main_seq,
	const sequence_vec & phylo_seqs,
	double threshold = BIOPSY_ANALYSE_THRESHOLD_DEFAULT );

/**
 * Calculate the maximal chain across the hit vectors.
 */
binding_hit::vec_ptr
analyse_max_chain(
	binding_hits_vec_ptr hit_array,
	unsigned max_box_limit );

/**
Find the pathway associated with a pssm
*/

std::string
get_pathway_for_pssm(
	const std::string & pssm_name);


/**
 * Build an SVG file.
 */
void
build_svg(
	const boost::filesystem::path & file,
	const std::string & title,
	const BIO_NS::seq_t & seq,
	double min_threshold,
	const binding_hit::vec & hits,
	size_t max_num_factors,
	bool show_labels,
	bool open,
	const binding_hit::vec * max_chain,
	const std::string & notes,
	double max_threshold = 0.0 );


unsigned
num_hits_at_threshold( const binding_hit::vec & hits, double threshold );


unsigned
num_boxes_for_threshold(
    binding_hits_vec_ptr hit_array,
    double threshold
);

typedef std::pair< double, unsigned > size_for_threshold;
typedef std::pair< size_for_threshold, size_for_threshold > threshold_pair;


threshold_pair
calculate_best_threshold(
    binding_hits_vec_ptr hit_array,
    unsigned num_boxes_limit
);

///**
// * We are given a range that contains counts and need to reduce each
// * count so that the product is less than some maximum. We average the
// * amout to remove in log-space.
// */
//template< typename CountsRange >
//double
//calculate_average_count_reduction(
//    const CountsRange & counts,
//    unsigned max_product
//) {
//    // check counts make sense
//    const double log_max_product = std::log( max_product );
//    const double log_product = std::log( product( counts ) );
//    //BOOST_ASSERT( log_product >= log_max_product );
//
//    // how much do we need to remove from each count?
//    const double log_to_remove = log_product - log_max_product;
//    return log_to_remove / boost::size( counts );
//}


/** Calculate how many to remove given a count and a log-scale reduction. */
unsigned
how_many_to_remove( unsigned count, double log_to_remove_per_count );


} //namespace biopsy

#endif //BIOPSY_ANALYSE_H_
