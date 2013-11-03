/**
@file

Copyright John Reid 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013

*/
# ifdef _MSC_VER
# pragma warning(disable : 4503)
# endif //_MSC_VER

#include "biopsy/defs.h"
#include "biopsy/analyse.h"
#include "biopsy/pssm.h"
#include "biopsy/sequence.h"
#include "biopsy/bifa.h"

#include "bio/bifa_analysis.h"
#include "bio/adjust_hits.h"
#include "bio/biobase_score.h"
#include "bio/biobase_filter.h"
#include "bio/binding_model.h"
#include "bio/biobase_binding_model.h"
#include "bio/pathway_associations.h"

USING_BIO_NS


#include "gsl/gsl_math.h"


namespace biopsy {
namespace detail {

binding_hit convert( const BIO_NS::BindingHit< BindingModel > & hit )
{
    return
        binding_hit(
            hit.get_binder()->get_name(),
            binding_hit_location(
                hit.get_position(),
                hit.get_length(),
                ! hit.is_complementary()),
            hit.get_p_binding() );
}

binding_hit::vec_ptr convert( const bifa_hits_t & bifa_hits )
{
    binding_hit::vec_ptr result( new binding_hit::vec );
    BOOST_FOREACH( const BIO_NS::BindingHit< BindingModel > & hit, bifa_hits )
    {
        result->push_back( convert( hit ) );
    }
    return result;
}

} //namespace detail






/**
 * Evaluate all the words in the sequence. Returns probability of binding at least once to sequence.
 */
template< typename Evaluator >
double
evaluate_words_in_sequence(
    const pssm_info &      info,
    const std::string &    pssm_name,
    const sequence &       seq,
    double                 threshold,
    const Evaluator &      evaluator,
    binding_hit::vec_ptr   result
) {
    double p_does_not_bind_anywhere = 1.0;
    size_t position = 0;
    const size_t size = info._pssm->size();
    for( sequence::const_iterator s = seq.begin(); true; ++s, ++position )
    {
        //is there enough sequence left to score this pssm?
        if( seq.end() - s < int( size ) )
        {
            break; //no
        }

#if ! defined( NDEBUG )
        const std::string dbg_word( s, s + size );
#endif // ! defined( NDEBUG )

        //are there any 'n's before the end of the pssm?
        sequence::const_iterator pssm_end = s + size;
        if( pssm_end != std::find_if( s, pssm_end, is_unknown_nucleotide() ) )
        {
            continue;
        }

        for( int i = 0; 2 != i; ++i ) // i=0 for positive strand, i=1 for negative strand
        {
            const bool is_positive_strand = (0 == i);
            const double p_binding = evaluator( s, is_positive_strand, position );
            if( p_binding >= threshold )
            {
                result->push_back(
                    binding_hit(
                        pssm_name,
                        binding_hit_location(
                            position,
                            size,
                            is_positive_strand ),
                        p_binding ) );

                p_does_not_bind_anywhere *= ( 1.0 - p_binding );
            }
        }
    }

    const double p_binds_somewhere = 1.0 - p_does_not_bind_anywhere;
    return p_binds_somewhere;
}


/// Functor that evaluates a word using the older score method.
struct evaluate_word_using_score {

    const pssm_info & info;

    evaluate_word_using_score( const pssm_info & info ) : info( info ) { }

    // Evaluate the word.
    double
    operator()(
        sequence::const_iterator s,
        bool is_positive_strand,
        size_t position
    ) const {
        return is_positive_strand
            ? get_p_binding_on_sequence( info, s )
            : get_p_binding_on_reverse_complement( info, s );
    }
};




/// Functor that evaluates a word using the BiFa method.
template< typename BgLikelihoods >
struct evaluate_word_using_bifa {

    typedef pssm_info::matrix_t                                             pssm_t;

    const BgLikelihoods &                                 bg_likelihoods;
    const pssm_t &                                          pssm_log_likelihoods;
    typename bifa::PssmTraits< pssm_t >::reverse_complement pssm_rev_comp_log_likelihoods;

    evaluate_word_using_bifa(
        const pssm_info & info,
        const BgLikelihoods & bg_likelihoods
    )
    : bg_likelihoods( bg_likelihoods )
    , pssm_log_likelihoods( info.get_log_likelihoods() )
    , pssm_rev_comp_log_likelihoods( bifa::pssm_reverse_complement( const_cast< pssm_t & >( pssm_log_likelihoods ) ) )
    { }

    // Evaluate the word.
    double
    operator()(
        sequence::const_iterator s,
        bool is_positive_strand,
        size_t position
    ) const {
        const double pssm_log_likelihood =
            is_positive_strand
                ? bifa::score_word< bifa::DnaAlphabet >(
                    pssm_log_likelihoods,
                    boost::make_transform_iterator( s, biopsy::bifa::convert_char_base_to_int() )
                )
                : bifa::score_word< bifa::DnaAlphabet >(
                    pssm_rev_comp_log_likelihoods,
                    boost::make_transform_iterator( s, biopsy::bifa::convert_char_base_to_int() )
                )
            ;
        const double bg_log_likelihood = bg_likelihoods.get_word_log_likelihood(
            position,
            boost::size( pssm_log_likelihoods )
        );

        // calculate the probability of binding using the odds ratio.
        return get_p_binding(
            pssm_parameters::singleton().binding_background_odds_prior
            * std::exp( pssm_log_likelihood - bg_log_likelihood )
        );
    }
};




double
score_pssm_on_sequence(
    const std::string & pssm_name,
    const sequence & seq,
    double threshold,
    binding_hit::vec_ptr result
) {
    const pssm_info & info = get_pssm( pssm_name );
    const pssm_parameters & params = pssm_parameters::singleton();

    return
        params.use_score
            ? evaluate_words_in_sequence(
                info,
                pssm_name,
                seq,
                threshold,
                evaluate_word_using_score( info ),
                result
            )
            : evaluate_words_in_sequence(
                    info,
                    pssm_name,
                    seq,
                    threshold,
                    evaluate_word_using_bifa<
                        bifa::uniform_sequence_likelihoods
                    >( info, bifa::uniform_sequence_likelihoods() ),
                    result
            )
        ;
}



void
biobase_score_pssm_on_sequence(
    const std::string & pssm_name,
    const sequence & seq,
    double threshold,
    binding_hit::vec_ptr result )
{
    const pssm_info info = get_pssm( pssm_name );

    unsigned position = 0;
    for( sequence::const_iterator s = seq.begin(); true; ++s, ++position )
    {
        //is there enough sequence left to score this pssm?
        if( seq.end() - s < int( info._pssm->size() ) )
        {
            break; //no
        }

        //are there any 'n's before the end of the pssm?
        sequence::const_iterator pssm_end = s + info._pssm->size();
        if( pssm_end != std::find_if( s, pssm_end, is_unknown_nucleotide() ) )
        {
            continue;
        }

        for( int i = 0; 2 != i; ++i )
        {
            const double biobase_score =
                0 == i
                    ? score( *(info._pssm), s )
                    : score_complement( *(info._pssm), s );
            const bool is_positive_strand = (0 == i);
            if( biobase_score >= threshold )
            {
                result->push_back(
                    binding_hit(
                        pssm_name,
                        binding_hit_location(
                            position,
                            info._pssm->size(),
                            is_positive_strand ),
                        biobase_score ) );
            }
        }
    }
}


binding_hit::vec_ptr
score_pssms_on_sequence(
    const string_vec_ptr & pssm_names,
    const sequence & seq,
    double threshold )
{
    binding_hit::vec_ptr result( new binding_hit::vec );
    BOOST_FOREACH( const std::string & pssm_name, *pssm_names )
    {
        score_pssm_on_sequence(
            pssm_name,
            seq,
            threshold,
            result );
    }

    return result;
}


binding_hit::vec_ptr
biobase_score_pssms_on_sequence(
    const string_vec_ptr & pssm_names,
    const sequence & seq,
    double threshold )
{
    binding_hit::vec_ptr result( new binding_hit::vec );
    BOOST_FOREACH( const std::string & pssm_name, *pssm_names )
    {
        biobase_score_pssm_on_sequence(
            pssm_name,
            seq,
            threshold,
            result );
    }

    return result;
}


/**
 * Abstract base class for phylogenetic adjusters.
 */
struct phylogenetic_adjuster {
    virtual ~phylogenetic_adjuster();

    /**
     * Accept the probability that the binder binds at least once to one of the
     * phylogenetic sequences.
     */
    virtual
    void
    accept_prob_phylo_binding( double p_binding ) = 0;


    /**
     * Adjust the strength of the hit in the central sequence given the previously accepted
     * probabilities of binding in the phylogenetic sequences.
     */
    virtual
    double
    adjust_hit_probability( unsigned num_sequences, double p_central ) = 0;
};


phylogenetic_adjuster::~phylogenetic_adjuster() { }


/**
 * Adjusts the probability of hits in the central sequence by averaging with the
 * probability that they bind anywhere in the related phylogenetic sequences.
 */
struct phylogenetic_adjuster_probability_averager
: phylogenetic_adjuster
{
    double log_sum;

    phylogenetic_adjuster_probability_averager() : log_sum( 0. ) { }
    virtual ~phylogenetic_adjuster_probability_averager() { }

    /**
     * Accept the probability that the binder binds at least once to one of the
     * phylogenetic sequences.
     */
    virtual
    void
    accept_prob_phylo_binding( double p_binding ) {
        log_sum += std::log( p_binding );
    }


    /**
     * Adjust the strength of the hit in the central sequence given the previously accepted
     * probabilities of binding in the phylogenetic sequences.
     */
    virtual
    double
    adjust_hit_probability( unsigned num_sequences, double p_central ) {
        //
        // Return geometric average of probabilities.
        //
        return
            std::exp(
                ( std::log( p_central ) + log_sum )
                / num_sequences
            );
    }
};


/**
 * Adjusts the probability of hits in the central sequence by averaging the weight
 * of evidence (Bayes factors) with the Bayes factors of the events that they bind
 * anywhere in the related phylogenetic sequences.
 */
struct phylogenetic_adjuster_bayes_averager
: phylogenetic_adjuster
{
    double log_sum; ///< The sum of the Bayes factors.
    const double prior_log_odds; ///< The prior log-odds of a binding site.
    const double min_log_bayes_factor; ///< The minimum log-Bayes factor that we will use. Designed to avoid problems with TFBSs missing in related sequences.

    phylogenetic_adjuster_bayes_averager(
        double log_prior_odds,
        double central_binding_p // The probability of binding in the central sequence
    )
    : log_sum( 0. )
    , prior_log_odds( log_prior_odds )
    , min_log_bayes_factor( calculate_min_log_bayes_factor( log_prior_odds, central_binding_p ) )
    { }
    virtual ~phylogenetic_adjuster_bayes_averager() { }

    /** Get the minimum log-Bayes factor we will use for the evidence from
     * related sequences. This is specified as a fraction of the evidence from
     * the central sequence.
     */
    static
    double
    calculate_min_log_bayes_factor( double log_prior_odds, double central_binding_p ) {
        const pssm_parameters & params = pssm_parameters::singleton();
        BOOST_ASSERT( 0. <= params.min_related_evidence_fraction );
        BOOST_ASSERT( params.min_related_evidence_fraction <= 1. );
        // if fraction is turned off, then the minimum evidence (log Bayes factor) does not apply
        if( ! params.min_related_evidence_fraction ) {
            return -std::numeric_limits< double >::max();
        } else {
            const double central_log_bayes_factor = prob_to_log_odds( central_binding_p ) - log_prior_odds;
            // if negative evidence in central sequence do not reduce it
            if( central_log_bayes_factor < 0. ) {
                return central_log_bayes_factor;
            } else {
                return central_log_bayes_factor * params.min_related_evidence_fraction;
            }
        }
    }

    /**
     * The log-odds given a probability.
     */
    static
    double
    prob_to_log_odds( double p ) {
        return std::log( p / ( 1. - p ) );
    }

    /**
     * The probability given the log-odds.
     */
    static
    double
    log_odds_to_prob( double log_odds ) {
        const double odds = std::exp( log_odds );
        BOOST_ASSERT( ! BIO_ISNAN( odds ) );
        if( BIO_FINITE( odds ) ) {
            return odds / ( 1. + odds );
        } else {
            return 1.;
        }
    }

    /**
     * Accept the probability that the binder binds at least once to one of the
     * phylogenetic sequences. Here we make an assumption that the length
     * of the phylogenetic sequence is smaller than 1/prior_odds
     */
    virtual
    void
    accept_prob_phylo_binding( double p_binding ) {
        BOOST_ASSERT( ! BIO_ISNAN( p_binding ) );
        const double posterior_log_odds = prob_to_log_odds( p_binding );
        const double log_bayes_factor = posterior_log_odds - prior_log_odds;
        log_sum += std::max( min_log_bayes_factor, log_bayes_factor );
        // BOOST_ASSERT( BIO_FINITE( log_sum ) );  allow infinite log sums
    }


    /**
     * Adjust the strength of the hit in the central sequence given the previously accepted
     * probabilities of binding in the phylogenetic sequences.
     */
    virtual
    double
    adjust_hit_probability( unsigned num_sequences, double p_central ) {
        const double central_log_odds = prob_to_log_odds( p_central );
        const double central_log_bayes_factor = central_log_odds - prior_log_odds;
        const double avg_log_bayes_factor = ( central_log_bayes_factor + log_sum ) / num_sequences;
        const double p_adjusted = log_odds_to_prob( prior_log_odds + avg_log_bayes_factor );
        BOOST_ASSERT( ! BIO_ISNAN( p_adjusted ) );
        return p_adjusted;
    }
};


phylo_sequences_result
score_pssms_on_phylo_sequences(
    string_vec_ptr pssm_names_arg,
    sequence_vec_ptr sequences,
    double threshold,
    double phylo_threshold,
    bool calculate_maximal_chain
) {
    //
    // We need at least one sequence
    //
    if( sequences->empty() ) {
        throw std::invalid_argument( "Need at least one sequence to score" );
    }

    //
    // the prior odds of a binding site
    //
    const pssm_parameters & params = pssm_parameters::singleton();
    const double prior_log_odds = std::log( params.binding_background_odds_prior );

    //
    // make a copy of the pssm names argument
    //
    string_vec_ptr pssm_names( new string_vec( *pssm_names_arg ) );

    //
    // A map from binders to phylogenetic adjusters
    //
    std::map< std::string, boost::shared_ptr< phylogenetic_adjuster > > phylo_adjusters;

    //
    // for each sequence score the pssms we are interested in
    //
    binding_hits_vec_ptr hit_array( new binding_hits_vec );
    bool is_first_sequence = true;
    BOOST_FOREACH( const sequence & s, *sequences ) {

        //push_back a hit vector for this sequence
        hit_array->push_back( binding_hit::vec_ptr( new binding_hit::vec ) );
        binding_hit::vec_ptr hits = hit_array->back();
        BOOST_FOREACH( const std::string & pssm_name, *pssm_names ) {

            try {
                const double binding_p =
                    score_pssm_on_sequence(
                        pssm_name,
                        s,
                        is_first_sequence
                            ? threshold
                            : phylo_threshold,
                        hits
                    );

                //
                // Phylogenetic adjustment stuff
                //
                if( is_first_sequence ) {
                    //
                    // Create a phylogenetic adjuster for this PSSM
                    //
                    if( params.avg_phylo_bayes ) {
                        phylo_adjusters[ pssm_name ].reset( new phylogenetic_adjuster_bayes_averager( prior_log_odds, binding_p ) );
                    } else {
                        phylo_adjusters[ pssm_name ].reset( new phylogenetic_adjuster_probability_averager );
                    }
                } else {
                    //
                    // Pass the probability of binding to the phylogenetic adjuster
                    //
                    phylo_adjusters[ pssm_name ]->accept_prob_phylo_binding( binding_p );
                }
            } catch( std::exception const & e ) {
                throw std::logic_error(
                    BIOPSY_MAKE_STRING(
                        "Problem scoring PSSM: "<<pssm_name<<": "<<e.what() ) );
            } catch( ... ) {
                throw std::logic_error(
                    BIOPSY_MAKE_STRING(
                        "Unknown problem scoring PSSM: "<<pssm_name ) );
            }
        }

        //
        // Rebuild set of pssms we are interested in. We will
        // not bother scoring PSSMs that we do not have hits for
        // in every sequence so far
        //
        pssm_names = get_binder_names( hits );

        // It won't be the first sequence next time around
        is_first_sequence = false;
    }

    //remove hits for pssms that are not in all sequences - pssm_names holds those pssms that are
    BOOST_FOREACH( binding_hit::vec_ptr & hits, *hit_array ) {
        //create new vector to hold filtered hits
        binding_hit::vec_ptr filtered_hits( new binding_hit::vec );

        //copy those hits that are in the pssm_names container
        BOOST_FOREACH( const binding_hit & hit, *hits ) {
            if( pssm_names->end() != std::find( pssm_names->begin(), pssm_names->end(), hit._binder_name ) ) {
                filtered_hits->push_back( hit );
            }
        }

        //replace original hits with filtered
        hits.swap( filtered_hits );
    }

    //calculate the maximal chain if we can and want to
    binding_hit::vec_ptr mc;
    if( calculate_maximal_chain ) {
        mc = analyse_max_chain(
            hit_array,
            pssm_parameters::singleton().max_chain_num_boxes_limit
        );
    }

    //
    // Adjust the hits for phylogenetic conservation
    //
    if( ! hit_array->empty() ) {
        BOOST_FOREACH( binding_hit & hit, *hit_array->front() ) {
            // adjust by estimate that binds in the phylo sequences
            BOOST_ASSERT( ! BIO_ISNAN( hit._p_binding ) );
            hit._p_binding = phylo_adjusters[ hit._binder_name ]->
                adjust_hit_probability( hit_array->size(), hit._p_binding );
            BOOST_ASSERT( ! BIO_ISNAN( hit._p_binding ) );
        }
    }

    phylo_sequences_result result( hit_array->front(), mc, hit_array );
    return result;
}

binding_hit::vec_ptr
analyse(
    const sequence & seq,
    double threshold )
{
    bifa_hits_t bifa_hits;
    {
        BiobasePssmFilter filter;
        score_all_biobase_pssms(
            make_sequence_scorer(
                seq.begin(),
                seq.end(),
                threshold,
                std::inserter( bifa_hits, bifa_hits.begin() )
            ),
            filter,
            Link2BiobaseBindingModel()
        );
    }

    return detail::convert( bifa_hits );
}



binding_hit::vec_ptr
analyse_phylo(
    const sequence & main_seq,
    const sequence_vec & phylo_seqs,
    double threshold )
{
    bifa_hits_t bifa_hits;
    {
        BiobasePssmFilter filter;
        score_all_biobase_pssms(
            make_sequence_scorer(
                main_seq.begin(),
                main_seq.end(),
                threshold,
                std::inserter( bifa_hits, bifa_hits.begin() )
            ),
            filter,
            Link2BiobaseBindingModel()
        );

        //raise the threshold to the power of the number of sequences
        const BIO_NS::float_t phylo_threshold = BIO_NS::float_t( gsl_pow_int( threshold, phylo_seqs.size() + 1 ) );

        adjust_hits(
            bifa_hits,
            phylo_seqs,
            phylo_threshold);
    }

    return detail::convert( bifa_hits );

}



std::string
get_pathway_for_pssm( const std::string & pssm_name ) {

/*    binding_hit::vec empty_chain;
    BiFaDetails details( hits, empty_chain );
    set_pathways(details);*/
    return "";

}


} //namespace biopsy

