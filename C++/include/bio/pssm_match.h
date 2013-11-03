#ifndef BIO_PSSM_MATCH_H_
#define BIO_PSSM_MATCH_H_


#include "bio/defs.h"
#include "bio/biobase_match.h"
#include "bio/environment.h"
#include "bio/biobase_likelihoods.h"
#include "bio/matrix_match.h"
#include "bio/score_sequence.h"


#include <boost/shared_ptr.hpp>

BIO_NS_START

/** Maps pssms to results/scores. */
template <class ResultCalculator>
struct PssmScoreMap : std::map<TableLink, hit_vec_t>
{
	PssmScoreMap(const ResultCalculator & result_calculator)
		: _result_calculator(result_calculator)
	{
	}

	template <class Pssm>
	const hit_vec_t &
	get_result(const Pssm & pssm)
	{
		//find the scores in the map
		iterator i = find(pssm.get_link());

		//did we find them?
		if (end() == i)
		{
			//no - create them

			//first insert empty result vector
			i = insert(make_pair(pssm.get_link(), hit_vec_t())).first;

			_result_calculator(pssm, inserter(i->second, i->second.begin()));
		}

		return i->second;
	}

	size_t
	get_num_hits()
	{
		size_t total_hits = 0;
		for (const_iterator j = begin(); end() != j; ++j)
		{
			total_hits += j->second.size();
		}
		return total_hits;
	}

	/** Ranks one particular hit against all others in the map. */
	size_t
	rank_hit(
		const TableLink & link, 
		float_t score) const
	{
		//find the hit vector for this link
		const_iterator i = find(link);
		if (end() == i)
		{
			throw std::logic_error( "Could not find data for link" );
		}
		const hit_vec_t & hits = i->second;

		size_t rank = 0; // the rank of our best_score over all the hits.

		//find the rank out of all the hits in the map
		for (const_iterator j = begin(); end() != j; ++j)
		{
			rank += rank_score(j->second, score);
		}

		return rank;
	}

	/** Inserts the pssms that reached the threshold into the insert iterator. */
	template <class InsIt>
	void
	get_pssms_over_threshold(
		float_t threshold,
		InsIt insert_iterator) const
	{
		for (const_iterator i = begin();
			end() != i;
			++i)
		{
			hit_vec_t::const_iterator h = i->second.begin();
			while (i->second.end() != h)
			{
				if (h->score > threshold)
				{
					break;
				}
				++h;
			}
			if (i->second.end() != h)
			{
				*insert_iterator++ = i->first;
			}
		}
	}
	

protected:
	ResultCalculator _result_calculator;
};





/** Calculates biobase scores over a sequence. */
struct BiobaseResultCalculator
{
	BiobaseResultCalculator(seq_t::const_iterator match_seq_begin, seq_t::const_iterator match_seq_end)
		: match_seq_begin(match_seq_begin)
		, match_seq_end(match_seq_end)
	{
	}

	template <class Pssm, class InsIt>
	void operator()(const Pssm & pssm, InsIt insert_it) const
	{
		//get all the scores for this pssm
		score_sequence(
			make_pssm(&pssm),
			match_seq_begin,
			match_seq_end,
			true,
			insert_it);
	}

protected:
	seq_t::const_iterator  match_seq_begin;
	seq_t::const_iterator  match_seq_end;
};
typedef PssmScoreMap <BiobaseResultCalculator> BiobaseScoresMap;


/** Calculates ott normalised scores over a sequence. */
struct BiobaseFilterResultCalculator
{
	BiobaseFilterResultCalculator(BiobaseScoresMap & biobase_scores_map)
		: biobase_scores_map (biobase_scores_map)
	{
	}

	template <class Pssm, class InsIt>
	void operator()(const Pssm & pssm, InsIt insert_it) const
	{
#ifdef _MSC_VER
		//find the threshold for this pssm
		MatrixMatch::map_t::const_iterator mm = MatrixMatchMap::singleton().find(pssm.get_link());
		const float_t threshold =
			MatrixMatchMap::singleton().end() == mm
				? 0.0f // cannot find threshold - use 0
				: mm->second.min_fp_threshold;
#else //_MSC_VER
		// It looks like MSVC never bothered to instantiate this function.
		// gcc tries to compile it but can't because MatrixMatchMap is not defined anywhere
		// We throw an error to at least catch this at runtime.
		const float_t threshold = 0.0;
		throw std::logic_error( "In uncompilable code!" );
#endif

		//find the biobase scores
		const hit_vec_t & biobase_scores = biobase_scores_map.get_result(pssm);

		//insert only the hits above the threshold
		for (hit_vec_t::const_iterator i = biobase_scores.begin();
			biobase_scores.end() != i;
			++i)
		{
			if (threshold < i->score)
			{
				*insert_it++ =
					Hit(
						i->score,
						i->position,
						i->complement);
			}
		}
	}

protected:
	BiobaseScoresMap & biobase_scores_map;
};
typedef PssmScoreMap <BiobaseFilterResultCalculator> BiobaseFilterScoresMap;


/** Calculates ott normalised scores over a sequence. */
struct OttNormalisedResultCalculator
{
	OttNormalisedResultCalculator(BiobaseScoresMap & biobase_scores_map)
		: biobase_scores_map (biobase_scores_map)
	{
	}

	template < class Pssm, class InsIt >
	void operator()(const Pssm & pssm, InsIt insert_it) const
	{
		//get the normaliser
		const BiobaseLikelihoods * normaliser =
			LikelihoodsCache::singleton().get_score_likelihoods( pssm.get_link(), true, true );
		if (0 == normaliser)
		{
			throw std::logic_error( BIO_MAKE_STRING( ("No normaliser for pssm ") + pssm.get_name() ) );
		}

		//find the biobase scores
		const hit_vec_t & biobase_scores = biobase_scores_map.get_result(pssm);

		//insert the normalised scores
		for (hit_vec_t::const_iterator i = biobase_scores.begin();
			biobase_scores.end() != i;
			++i)
		{
			BOOST_ASSERT( 0 <= i->score );
			BOOST_ASSERT( i->score <= 1.0 );
			*insert_it++ =
				Hit(
					1.0f - get_likelihood( *normaliser, i->score ),
					i->position,
					i->complement);
		}
	}

protected:
	BiobaseScoresMap & biobase_scores_map;
};
typedef PssmScoreMap <OttNormalisedResultCalculator> OttNormalisedScoresMap;


/** Calculates Bayesian scores over a sequence. */
template <bool or_better>
struct BayesianResultCalculator
{
	BayesianResultCalculator(BiobaseScoresMap & biobase_scores_map)
		: biobase_scores_map (biobase_scores_map)
	{
	}

	template <class Pssm, class InsIt>
	void operator()(const Pssm & pssm, InsIt insert_it) const
	{
		using namespace boost::numeric;

		//define the priors - H0 is non-binding hypothesis, H1 is binding hypothesis
		//if in 400 conserved bases (and complementary) we expect to see 4 binding sites after testing 1000 pssms
		//our prior should be:
		static const float_t p_H1 = BioEnvironment::singleton().get_tf_binding_prior();
		static const float_t p_H0 = float_t(1.0) - p_H1;
		BOOST_ASSERT( float_t(0.0) < p_H1 && p_H1 < float_t(1.0) );

		//get the normaliser
		const BiobaseLikelihoods * background =
			LikelihoodsCache::singleton().get_score_likelihoods(pssm.get_link(), true, or_better);
		const BiobaseLikelihoods * binding =
			LikelihoodsCache::singleton().get_score_likelihoods(pssm.get_link(), false, or_better);
		if (0 == background || 0 == binding)
		{
			throw BIO_MAKE_STRING("No normaliser for pssm: " << pssm.get_name());
		}

		//find the biobase scores
		const hit_vec_t & biobase_scores = biobase_scores_map.get_result(pssm);

		//insert the normalised scores
		for (hit_vec_t::const_iterator i = biobase_scores.begin();
			biobase_scores.end() != i;
			++i)
		{
			const float_t p_D_given_H1 = get_likelihood(*binding, i->score);
			BOOST_ASSERT( in( p_D_given_H1, interval< float_t >( 0.0f, 1.0f ) ) );

			//this could be 0
			const float_t p_D_given_H0 = get_likelihood(*background, i->score);
			BOOST_ASSERT( in( p_D_given_H0, interval< float_t >( 0.0f, 1.0f ) ) );
			
			//so this could be infinite
			const float_t bayes_factor = (p_D_given_H1 * p_H1) / (p_D_given_H0 * p_H0);

			//in which case this should be 1
			const float_t p_H1_given_D =
				BIO_FINITE( bayes_factor )
					? bayes_factor / (1.0f + bayes_factor)
					: 1.0f;

			//our prob should be between 0 and 1
			if ( 0.0f > p_H1_given_D || p_H1_given_D > 1.0f )
			{
				throw std::logic_error( "Score is out of range" );
			}
			BOOST_ASSERT( in( p_H1_given_D, interval< float_t >( 0.0f, 1.0f ) ) );

			//put our hit in the results
			*insert_it++ =
				Hit(
					p_H1_given_D,
					i->position,
					i->complement);
		}
	}

protected:
	BiobaseScoresMap & biobase_scores_map;
};
typedef PssmScoreMap <BayesianResultCalculator<false> > BayesianScoresMap;
typedef PssmScoreMap <BayesianResultCalculator<true> > BayesianOrBetterScoresMap;


struct Scorers
{
	typedef boost::shared_ptr<Scorers> ptr_t;

	Scorers(seq_t::const_iterator begin, seq_t::const_iterator end)
		: biobase_calc(begin, end)
		, biobase_scores(biobase_calc)
		, ott_calc(biobase_scores)
		, ott_normalised_scores(ott_calc)
		, bayesian_calc(biobase_scores)
		, bayesian_scores(bayesian_calc)
		, bayesian_or_better_calc(biobase_scores)
		, bayesian_or_better_scores(bayesian_or_better_calc)
		, biobase_filter_calc(biobase_scores)
		, biobase_filter_scores(biobase_filter_calc)
	{
	}

	BiobaseResultCalculator biobase_calc;
	BiobaseScoresMap biobase_scores;

	OttNormalisedResultCalculator ott_calc;
	OttNormalisedScoresMap ott_normalised_scores;

	BayesianResultCalculator<false> bayesian_calc;
	BayesianScoresMap bayesian_scores;

	BayesianResultCalculator<true> bayesian_or_better_calc;
	BayesianOrBetterScoresMap bayesian_or_better_scores;

	BiobaseFilterResultCalculator biobase_filter_calc;
	BiobaseFilterScoresMap biobase_filter_scores;
};

BIO_NS_END



#endif //BIO_PSSM_MATCH_H_
