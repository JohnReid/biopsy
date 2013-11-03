/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/pssm_likelihood.h"


#include <boost/numeric/interval.hpp>


#include <iostream>
#include <sstream>




BIO_NS_START
/** Calculates the likelihoods of certain biobase scores given sequences that match the pssm. */ 
void
pssm_likelihood(
    const Pssm & pssm,
    const size_t max_map_size,
    BiobaseLikelihoods & result,
    float_t pseudo_count,
    bool verbose)
{
    using namespace boost::numeric;


    typedef double prob_t; //using float loses accuracy

    /** Maps scores to probabilities. */
    typedef std::map<prob_t, prob_t> probs_t;
    boost::shared_ptr<probs_t> probs1(new probs_t);
    boost::shared_ptr<probs_t> probs2(new probs_t);
    probs_t * last_probs = probs1.get(); //these alternate
    probs_t * new_probs = probs2.get();
    (*last_probs)[0.0] = std::log(1.0); //initialise with score of 0 has likelihood 1.0 at the start of the algorithm

    //for each row
    size_t row_idx = 0;
    float_t min_score = 0.0;
    float_t max_score = 0.0;
    for (Pssm::const_iterator entry = pssm.begin();
        pssm.end() != entry;
        ++entry)
    {
        max_score += entry->get_max();
        min_score += entry->get_min();

        new_probs->clear(); //start afresh

        //go to the next nucleotide if no observations
        if (0 == entry->get_num_observations())
        {
            continue;
        }

        //if we have too many scores we need to quantise them
        while (last_probs->size() > max_map_size)
        {
            if (verbose)
            {
                std::cout << "Reducing map size from " <<  last_probs->size();
            }

            const prob_t min_score = last_probs->begin()->first;
            const prob_t max_score = last_probs->rbegin()->first;
            const prob_t quantum = (max_score - min_score) / max_map_size;

            //for each pair of consecutive scores
            for (probs_t::iterator p1 = last_probs->begin(); last_probs->end() != p1; )
            {
                //get the next score
                probs_t::iterator p2 = p1; ++p2;
                if (last_probs->end() == p2)
                {
                    break;
                }

                const prob_t score1 = p1->first;
                const prob_t score2 = p2->first;

                //are the two scores close?
                BOOST_ASSERT(score2 > score1);
                if (score2 - score1 < quantum)
                {
                    //we remove them both and insert a replacement for the average score that is as likely
                    //as both together
                    probs_t::iterator to_erase1 = p1;
                    probs_t::iterator to_erase2 = p2;

                    //insert the new value
                    const prob_t prob1 = std::exp(p1->second);
                    const prob_t prob2 = std::exp(p2->second);
                    const prob_t prob_sum = prob1 + prob2;
                    //const prob_t relative_prob1 = prob1 / prob_sum;
                    const prob_t relative_prob2 = prob2 / prob_sum;
                    const prob_t avg_score = score1 + (score2 - score1) * relative_prob2; //had to be careful here with floating point underflow
                    BOOST_ASSERT(score1 <= avg_score);
                    BOOST_ASSERT(avg_score <= score2);

                    //move on before we insert the new value
                    p1 = p2; ++p1;

                    //erase the old values
                    last_probs->erase(to_erase1);
                    last_probs->erase(to_erase2);

                    //is the average score in the map already?
                    probs_t::iterator already_in_map = last_probs->find( avg_score );
                    const prob_t p_already_seen = last_probs->end() == already_in_map ? 0.0 : std::exp( already_in_map->second );

                    //insert the new value
                    const prob_t log_prob_sum = std::log( prob_sum + p_already_seen );
                    //std::pair<probs_t::iterator, bool> insert_result =
                        last_probs->insert( probs_t::value_type( avg_score, log_prob_sum ) );
                }
                else
                {
                    //move on
                    p1 = p2; ++p1;
                }
            }

            if (verbose)
            {
                std::cout << " to " <<  last_probs->size() << std::endl;
            }

        } //while map is too large

        //for each nucleotide
        for (const char * n = "acgt"; 0 != *n; ++n)
        {
            //likelihood in the pssm for this nucleotide in this position? Use a psuedo count of 1
            const prob_t prob_this_nucleo = entry->get_freq( *n, pseudo_count );

            //if none, cannot contribute to the total probabilities
            if (0.0 == prob_this_nucleo)
            {
                continue;
            }
            //BOOST_ASSERT(0 != entry->get_count(*n)); - this was only true before pseudo counts used

            const prob_t log_prob_this_nucleo = std::log(prob_this_nucleo);
            const prob_t score_this_nucleo = entry->get_score(*n);

            //for each score we have already achieved
            for (probs_t::const_iterator p = last_probs->begin(); last_probs->end() != p; ++p)
            {
                const prob_t new_score = p->first + score_this_nucleo;
                const prob_t new_log_likelihood = p->second + log_prob_this_nucleo;

                //update the new probability map

                //first try to insert
                std::pair<probs_t::iterator, bool> insert_result =
                    new_probs->insert(probs_t::value_type(new_score, new_log_likelihood));
                
                //did it actually insert?
                if (! insert_result.second)
                {
                    //no - so update the existing element
                    insert_result.first->second =
                        std::log(std::exp(insert_result.first->second) + std::exp(new_log_likelihood));
                }
            }

        } //for each nucleotide

        if (verbose)
        {
            //for debugging
            prob_t total_prob = 0.0;
            for (probs_t::const_iterator p = new_probs->begin(); new_probs->end() != p; ++p)
            {
                total_prob += std::exp(p->second);
            }
            std::cout
                << "Row " << row_idx << " has prob map size " << new_probs->size()
                << " and total probability: " << total_prob
                << std::endl;

        }

        std::swap(new_probs, last_probs);

    } //for each row

    //initialise the result vector
    std::fill(result.begin(), result.end(), float_t(0.0));

    //for each score in the last_probs map, place in the likelihoods vector
    for (probs_t::const_iterator p = last_probs->begin(); last_probs->end() != p; ++p)
    {
        //normalise the score
        const prob_t score = p->first;
        prob_t normalised_score = (score - min_score) / (max_score - min_score);
        BOOST_ASSERT( in( normalised_score, interval< prob_t >( -0.0001, 1.0001 ) ) );
        normalised_score = std::max( prob_t( 0.0 ), normalised_score );
        normalised_score = std::min( prob_t( 1.0 ), normalised_score );

        //what is the probability
        const prob_t log_likelihood = p->second;
        const prob_t prob = std::exp(log_likelihood);

        result[ get_biobase_score_index( result.size(), float_t( normalised_score ) ) ] += float_t( prob );
    }

    //normalise result vector
    for (BiobaseLikelihoods::iterator i = result.begin();
        result.end() != i;
        ++i)
    {
        BOOST_ASSERT( in( *i, interval< float_t >( -0.0001f, 1.0001f ) ) );
        *i = std::max( float_t( 0.0 ), *i );
        *i = std::min( float_t( 1.0 ), *i );
    }
}

BIO_NS_END
