
#ifndef BIO_RUN_MATCH_H_
#define BIO_RUN_MATCH_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/match_hit.h"
#include "bio/sequence.h"
#include "bio/biobase_filter.h"
#include "bio/pssm_match.h"
#include "bio/phylogenetics.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"

#include <boost/filesystem/path.hpp>

#include <vector>
#include <set>

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>



namespace boost { namespace program_options {
    //forward decl
    class options_description;
} }



BIO_NS_START

enum ScoreAlgorithm
{
    OTT_SCORE_ALGORITHM,
    BAYESIAN_SCORE_ALGORITHM
};


struct MatchParams
{
    MatchParams(
        float_t t,
        ScoreAlgorithm alg = OTT_SCORE_ALGORITHM,
        bool or_better = false)
        : threshold( t )
        , score_algorithm( alg )
        , or_better( or_better )
    {
    }

    float_t threshold;
    ScoreAlgorithm score_algorithm;
    bool or_better;
};




struct MatchResults
{
    MatchResults();

    MatchResults(
        TableLink link,
        Hit result);

    TableLink link;
    Hit result;
    size_t number;

    bool operator<(const MatchResults & rhs) const
    {
        if (result == rhs.result)
        {
            if (link == rhs.link)
            {
                return number < rhs.number;
            }
            else
            {
                return link < rhs.link;
            }
        }
        else
        {
            return result < rhs.result;
        }
    }

protected:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & link;
        ar & result;
        ar & number;
    }
};
typedef std::vector<MatchResults> match_result_vec_t;
typedef boost::shared_ptr< match_result_vec_t > match_result_vec_ptr_t;
void sort_by_position(match_result_vec_t & matches);
std::ostream & operator<<(std::ostream & os, const MatchResults & match);


/** Compare based initially on position, then link. */
struct MatchResultsPositionLessThan
{
    bool
    operator()(const MatchResults & h1, const MatchResults & h2) const
    {
        return
            h1.result.position < h2.result.position
            ||
            (
                h1.result.position == h2.result.position
                &&
                h1.link < h2.link
            );
    }
};

inline
std::ostream &
operator<<(std::ostream & os, const MatchResults & details) {
    os
        << details.link << ","
        << details.result;
    return os;
}



/**
Score all the pssm on the sequence between match_seq_begin and match_seq_end and put the results in result_insert_it.
Returns an estimate that the pssm binds the sequence.
*/
template <class Pssm, class ResultInsIt>
float_t
score_pssm(
    const Pssm & pssm,
    seq_t::const_iterator match_seq_begin,
    seq_t::const_iterator match_seq_end,
    const MatchParams & params,
    ResultInsIt & result_insert_it)
{
    Scorers scorers(match_seq_begin, match_seq_end);

    //const hit_vec_t & hit_results = scorers.bayesian_scores.get_result(pssm);
    const hit_vec_t & hit_results =
        OTT_SCORE_ALGORITHM == params.score_algorithm
            ? scorers.ott_normalised_scores.get_result(pssm)
            : ( params.or_better
                ? scorers.bayesian_scores.get_result(pssm)
                : scorers.bayesian_or_better_scores.get_result(pssm));

    float_t p_does_not_bind = 1.0;
    for (hit_vec_t::const_iterator i = hit_results.begin();
        hit_results.end() != i;
        ++i)
    {
        if (i->score > params.threshold)
        {
            p_does_not_bind *= ( float_t( 1.0 ) - i->score );

            MatchResults result(pssm.get_link(), *i);
            assert(result.result.position < (int) (match_seq_end - match_seq_begin)); //make sure the position is not too high
            *result_insert_it++ = result;
        }
    }
    return float_t( 1.0 ) - p_does_not_bind;
}


/** Score all the pssms between begin and end and put the results in result_insert_it. */
template <class PssmIt, class ResultInsIt>
size_t
score_pssms(
    PssmIt pssm_begin,
    PssmIt pssm_end,
    seq_t::const_iterator match_seq_begin,
    seq_t::const_iterator match_seq_end,
    const MatchParams & params,
    ResultInsIt & result_insert_it)
{
    size_t num_pssms_matched = 0;
    for ( ; pssm_begin != pssm_end; ++pssm_begin)
    {
        try
        {
            score_pssm(
                *(pssm_begin->second.get()),
                match_seq_begin,
                match_seq_end,
                params,
                result_insert_it);

            ++num_pssms_matched;
        }
        catch (const std::exception & ex)
        {
            std::cerr << "Could not score " << pssm_begin->second->get_name() << ": " << ex.what() << std::endl;
        }
        catch (const std::string & ex)
        {
            std::cerr << "Could not score " << pssm_begin->second->get_name() << ": " << ex << std::endl;
        }
        catch (const char * ex)
        {
            std::cerr << "Could not score " << pssm_begin->second->get_name() << ": " << ex << std::endl;
        }
        catch (...)
        {
            std::cerr << "Could not score " << pssm_begin->second->get_name() << std::endl;
        }
    }
    return num_pssms_matched;
}

template <class SeqIt, class ResultInsIt, class ConsFilter, class MatrixFilter>
void
pssm_match(
    SeqIt seq_begin,
    SeqIt seq_end,
    float_t threshold,
    ScoreAlgorithm algorithm,
    bool use_or_better,
    ConsFilter consensus_filter,
    MatrixFilter matrix_filter,
    ResultInsIt result_inserter)
{
    typedef boost::filter_iterator<MatrixFilter, Matrix::map_t::const_iterator> matrix_filter_it;
    typedef boost::filter_iterator<ConsFilter, Site::map_t::const_iterator> site_filter_it;

    const matrix_filter_it matrices_begin(matrix_filter, BiobaseDb::singleton().get_matrices().begin(), BiobaseDb::singleton().get_matrices().end());
    const matrix_filter_it matrices_end(matrix_filter, BiobaseDb::singleton().get_matrices().end(), BiobaseDb::singleton().get_matrices().end());

    const site_filter_it sites_begin(consensus_filter, BiobaseDb::singleton().get_sites().begin(), BiobaseDb::singleton().get_sites().end());
    const site_filter_it sites_end(consensus_filter, BiobaseDb::singleton().get_sites().end(), BiobaseDb::singleton().get_sites().end());

    if (matrices_begin == matrices_end && sites_begin == sites_end)
    {
        throw std::logic_error( "No PSSMs match filter" );
    }

    MatchParams params(threshold, algorithm, use_or_better);
    {
        //std::cout << "Scoring matrices" << std::endl;
        //const size_t num_matched =
            score_pssms(
                matrices_begin,
                matrices_end,
                seq_begin,
                seq_end,
                params,
                result_inserter);
//                ostream_iterator<MatchResults>(cout, "\n"));
        //std::cout << "Scored " << num_matched << " matrices" << std::endl;
    }

    {
        //std::cout << "Scoring sites" << std::endl;
        //const size_t num_matched =
            score_pssms(
                sites_begin,
                sites_end,
                seq_begin,
                seq_end,
                params,
                result_inserter);
//                ostream_iterator<MatchResults>(cout, "\n"));
        //std::cout << "Scored " << num_matched << " sites" << std::endl;
    }
}



void
pssm_match(
    const seq_t & match_seq,
    float_t threshold,
    ScoreAlgorithm algorithm,
    const boost::filesystem::path & file,
    const std::string & title,
    bool show_labels);


template <typename HitIt>
float_t
estimate_binding_prob(
    HitIt hit_begin,
    HitIt hit_end)
{
    float_t product = 1.0;
    for (HitIt hit = hit_begin; hit_end != hit; ++hit)
    {
        product *= float_t( 1.0f - hit->result.score );
    }
    return float_t(1.0 - product);
}

template <typename Exponent>
double
power(double x, Exponent y)
{
    return
        0 == x
            ? 0
            : gsl_sf_exp(y * gsl_sf_log(x));
}

/** Adjusts hits for occurence in phylogenetically conserved sequence. */
template <typename ResultIt>
void
adjust_hits(
    ResultIt results_begin,
    ResultIt results_end,
    seq_t::const_iterator seq_begin,
    seq_t::const_iterator seq_end,
    float_t threshold,
    ScoreAlgorithm algorithm)
{
    const MatchParams params(threshold, algorithm);

    //maintain a set of already adjusted PSSMs
    std::set<TableLink> already_adjusted;

    //for each result
    for (ResultIt result = results_begin; results_end != result; ++result)
    {
        //check we haven't checked this PSSM before
        if (already_adjusted.end() != already_adjusted.find(result->link))
        {
            //we have so ignore it
            continue;
        }

        //we need to estimate likelihood that it binds in phylogenetic sequence so run PSSM over sequence.
        match_result_vec_t phylo_hits;
        std::insert_iterator<bio::match_result_vec_t> inserter(phylo_hits, phylo_hits.begin());
        switch (result->link.table_id)
        {
        case MATRIX_DATA:
            score_pssm(
                *BiobaseDb::singleton().get_entry<MATRIX_DATA>(result->link),
                seq_begin,
                seq_end,
                params,
                inserter);
            break;

        case SITE_DATA:
            score_pssm(
                *BiobaseDb::singleton().get_entry<SITE_DATA>(result->link),
                seq_begin,
                seq_end,
                params,
                inserter);
            break;

        default:
            throw std::logic_error( "Unknown pssm type" );
        }

        const float_t binding_prob = estimate_binding_prob(phylo_hits.begin(), phylo_hits.end());

        //adjust all the hits for this PSSM
        for (ResultIt to_adjust = result; results_end != to_adjust; ++to_adjust)
        {
            if (to_adjust->link == result->link)
            {
                //to_adjust->result.score *= power(binding_prob, 5.0 * phylo_seq.conservation * phylo_seq.conservation);
                to_adjust->result.score *= binding_prob;
            }
        }

        //add to the list of already adjusted PSSMs
        already_adjusted.insert(result->link);
    }
}

/** Adjusts hits for occurence in phylogenetically conserved sequence. */
template <
    typename ResultIt,
    typename SeqRange >
void
adjust_hits_for_phylo_sequences(
    ResultIt results_begin,
    ResultIt results_end,
    const SeqRange & sequences,
    float_t threshold,
    ScoreAlgorithm algorithm)
{
    unsigned num_seqs = 0;
    BOOST_FOREACH( const typename boost::range_value< SeqRange >::type & seq, sequences )
    {
        adjust_hits(
            results_begin,
            results_end,
            seq.begin(),
            seq.end(),
            threshold,
            algorithm );
        ++num_seqs;
    }
    std::cout << "Adjusted hits for " << num_seqs << " sequences" << std::endl;

    for( ; results_begin != results_end; ++results_begin )
    {
        //adjust hits - raise to 1/#seqs
        results_begin->result.score = 
            float_t(
                ( 0.0 == results_begin->result.score ) ? 0.0 : exp( log( results_begin->result.score ) / ( num_seqs + 1 ) ) );
    }
}

/** Adjust hits by power. */
struct RaiseHitToPower
{
    double power;

    RaiseHitToPower(double power);

    void operator()(MatchResults & hit) const;
};

/** Is a hit above a threshold? */
struct HitAboveThreshold
{
    float_t threshold;

    HitAboveThreshold(float_t threshold);

    bool operator()(const MatchResults & results) const;
};



BIO_NS_END


#endif //BIO_RUN_MATCH_H_
