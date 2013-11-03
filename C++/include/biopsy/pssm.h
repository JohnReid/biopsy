/**
@file

Copyright John Reid 2006, 2013

*/

#ifndef BIOPSY_PSSM_H_
#define BIOPSY_PSSM_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "biopsy/defs.h"

#include "bio/singleton.h"
#include "bio/pssm.h"


namespace biopsy {



/**
Defines a distribution/score for each base.
*/
class nucleo_dist
    : boost::equality_comparable< nucleo_dist >
{
public:
    typedef std::vector< nucleo_dist > vec;
    typedef boost::shared_ptr< vec > vec_ptr;

protected:
    /** A values for each nucleotide. */
    typedef boost::array< double, 4 > values;

    values _values;

public:
    nucleo_dist(
        double a = 0.0,
        double c = 0.0,
        double g = 0.0,
        double t = 0.0 );

    nucleo_dist(const nucleo_dist & counts,double pseudo_count);

    /** argument is 0 for a, 1 for c, ... */
    double get( unsigned nucleo ) const;
    void set( unsigned nucleo, double value );

    double get_total( ) const;

    /** argument is 0 for a, 1 for c, ... */
    double get_freq( unsigned nucleo ) const;

    bool operator==( const nucleo_dist & rhs ) const;

    friend class boost::serialization::access;
    template< typename  Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
        for( unsigned i = 0; _values.size() != i; ++i )
        {
            ar & _values[ i ];
        }
    }
};

nucleo_dist::vec operator + (const nucleo_dist::vec & v, double i);

nucleo_dist uniform_nucleo_dist();
nucleo_dist get_dist_for_iupac( char s );

/**
A pssm is just a score for each base at each position.
*/
typedef nucleo_dist::vec pssm;
typedef nucleo_dist::vec_ptr pssm_ptr;
typedef std::vector< pssm_ptr > pssm_ptr_vec;
typedef boost::shared_ptr< pssm_ptr_vec > pssm_ptr_vec_ptr;


/** Convert a biopsy pssm to a bio pssm. */
BIO_NS::Pssm make_transfac_pssm( const pssm & _pssm );

/** Get the IUPAC sequence for a PSSM. */
std::string calculate_iupac( const pssm & _pssm );

pssm_ptr
create_pssm(
    const nucleo_dist::vec & dists );

double
score_pssm(
    pssm_ptr pssm,
    const sequence & s );


double
score(
    const pssm & pssm,
    sequence::const_iterator s_begin );

double
score_complement(
    const pssm & pssm,
    sequence::const_iterator s_begin );


/**
A vector quantising the probabilities of any given score in [0,1].
*/
typedef std::vector< double > likelihoods;
typedef boost::shared_ptr< likelihoods > likelihoods_ptr;

/**
Get an index into a likelihoods vector given a score.
*/
unsigned get_likelihood_index( unsigned size, double score );


/**
Get the likelihoods of a score.
*/
double get_likelihood( likelihoods_ptr likelihoods, double score );


/**
Calculate the likelihoods of scores or greater.
*/
likelihoods_ptr accumulate_likelihoods( likelihoods_ptr likelihoods );


/**
Get the odds ratio for the binding and background hypotheses given a score.
*/
double
get_odds_ratio(
    double score,
    likelihoods_ptr _binding_dist,
    likelihoods_ptr _background_dist );


/**
Get the prob of binding given the odds ratio.
*/
double
get_p_binding(
    double odds_ratio );


/**
Get the odds ratio given the prob of binding.
*/
double
get_odds_ratio_from_p_binding(
    double p_binding );


/**
Normalises so that likelihoods add to 1. Returns value they used to sum to.
*/
double
normalise_likelihoods( likelihoods & result );


/**
Calculates the likelihoods of certain biobase scores given sequences that match the distribution.
*/
void
calculate_pssm_likelihoods(
    const pssm & pssm,
    const nucleo_dist::vec & distribution,
    likelihoods & result,
    const size_t max_map_size,
    bool verbose = false );




/**
A pssm together with likelihoods under its own distribution and a simple
background distribution.
*/
struct pssm_info
{
    //    Dists include the effect of pseudo counts.  _counts do not
    nucleo_dist::vec _counts;

    typedef boost::multi_array< double, 2 >        matrix_t;     ///< multi_array matrix type to store log likelihoods
    typedef boost::shared_ptr< matrix_t >          matrix_ptr;   ///< Pointer to matrix

    nucleo_dist::vec          _dists;
    double                    _pseudo_count;
    int                       _number_of_sites;
    pssm_ptr                  _pssm;
    likelihoods_ptr           _binding_dist;
    likelihoods_ptr           _background_dist;
    mutable likelihoods_ptr   _cumulative_binding_dist;
    mutable likelihoods_ptr   _cumulative_background_dist;
    mutable matrix_ptr        _log_likelihoods;

    pssm_info(
        const nucleo_dist::vec & counts = nucleo_dist::vec(),
        double pseudo_count = 0.,
        int number_of_sites = 1,
        pssm_ptr p = pssm_ptr(),
        likelihoods_ptr binding_dist = likelihoods_ptr(),
        likelihoods_ptr background_dist = likelihoods_ptr() );

    likelihoods_ptr get_dist( bool binding, bool cumulative ) const;

    /// Get the log likelihoods used in the BiFa algorithm
    const matrix_t & get_log_likelihoods() const;

    friend class boost::serialization::access;
    template< typename  Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & _counts;
        ar & _pseudo_count;
        ar & _number_of_sites;
        ar & _pssm;
        ar & _binding_dist;
        ar & _background_dist;
    }
};




/**
Creates/gets a PSSM (and its info) based on its id. This should match the following regex for
a transfac pssm... [MR][0-9][0-9][0-9][0-9][0-9]
*/
const pssm_info &
get_pssm( const std::string & id );

/** Get the PSSM's name based on its id. */
std::string get_pssm_name( const std::string & id );

/** Get the PSSM's URL based on its id. Returns "" if none. */
std::string get_pssm_url( const std::string & id );


/** Calculate likelihoods of scores under the PSSM's distribution. */
likelihoods_ptr
calculate_likelihoods_under_pssm( pssm_ptr _pssm, const nucleo_dist::vec & dists );


/** Calculate likelihoods of scores under a uniform background distribution. */
likelihoods_ptr
calculate_likelihoods_under_background( pssm_ptr _pssm );


/**
Calculate p(binding) given the pssm score using the p-value method.
*/
double
get_p_binding_using_p_value(
    double score,
    likelihoods_ptr likelihoods );



/**
Calculate p(binding) given the pssm score.
*/
double
get_p_binding_from_score(
    const pssm_info & p,
    double score );



/**
Scores a pssm on the sequence and adjusts for distributions.
*/
double
get_p_binding_on_sequence(
    const pssm_info & p,
    sequence::const_iterator s_begin );


/**
Scores a pssm on the reverse complement of the sequence and adjusts for distributions.
*/
double
get_p_binding_on_reverse_complement(
    const pssm_info & p,
    sequence::const_iterator s_begin );



/**
Add a pssm to the cache.
*/
bool
add_pssm_to_cache( const std::string & name, const pssm_info & pssm );


/**
Serialises the pssm cache.
*/
void
save_pssm_cache_state( );


/**
 * Clears the pssm cache. We may want to use this when testing parameter values.
 */
void
clear_pssm_cache( );


/**
 * Parameters for construction and scoring of pssms.
 */

struct pssm_parameters
    : BIO_NS::UserSingleton< pssm_parameters >
{
    ADD_STATIC_SINGLETON_VARIABLE( double,    pseudo_counts ) ///< The pseudo-counts to use when building PSSMs
    ADD_STATIC_SINGLETON_VARIABLE( unsigned,  likelihoods_size ) ///< The resolution of our distribution over PSSM scores.
    ADD_STATIC_SINGLETON_VARIABLE( unsigned,  calculate_likelihoods_map_size ) ///< The resolution of our map we use to calculate PSSM score likelihoods.
    ADD_STATIC_SINGLETON_VARIABLE( double,    binding_background_odds_prior ) ///< The prior odds of binding.
    ADD_STATIC_SINGLETON_VARIABLE( bool,      use_cumulative_dists ) ///< Use cumulative score distributions as opposed to exact.
    ADD_STATIC_SINGLETON_VARIABLE( bool,      use_p_value ) ///< Use a p-value based scoring method rather than comparing binding and background score likelihoods.
    ADD_STATIC_SINGLETON_VARIABLE( bool,      use_score ) ///< Use a PSSM scoring based scheme rather than the BiFA log-likelihood ratio algorithm
    ADD_STATIC_SINGLETON_VARIABLE( bool,      avg_phylo_bayes ) ///< Average Bayes factors instead of probabilities when adjusting for phylogenetic sequences.
    ADD_STATIC_SINGLETON_VARIABLE( unsigned, max_chain_num_boxes_limit ) ///< Limit on number of boxes used to calculate the maximal chain.
    ADD_STATIC_SINGLETON_VARIABLE( double,    min_related_evidence_fraction ) ///< The fraction of the evidence for binding from the central sequence that acts as a minimum for the related sequences.

    pssm_parameters();
};



} //namespace biopsy

BOOST_CLASS_TRACKING( biopsy::likelihoods, boost::serialization::track_always )

#endif //BIOPSY_PSSM_H_

