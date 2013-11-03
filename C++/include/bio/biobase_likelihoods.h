#ifndef BIO_BIOBASE_LIKELIHOOD_H_
#define BIO_BIOBASE_LIKELIHOOD_H_

#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/biobase_match.h"
#include "bio/singleton.h"

#include <boost/shared_ptr.hpp>

#include <string>
#include <vector>
#include <map>






BIO_NS_START



//forward decl
struct BiobaseDb;



/** Gets the index of a given score in [0,1]. */
size_t get_biobase_score_index(size_t size, float_t score);

/**
 * Normalises biobase match scores. Splits the [0,1] range up into quanta and stores a likelihood of the
 * biobase range falling in each quantum.
 */
typedef std::vector<float_t> BiobaseLikelihoods;

/** Quantises the biobase scores and counts the number of hits in each quantum. In some contexts
the entries represent cumulative counts. */
typedef std::vector<size_t> BiobaseCounts;


/** Turn the counts of quantised scores into cumulative counts. */
BiobaseCounts
create_cumulative_from_counts(const BiobaseCounts & counts);


/** Given a vector of counts, calculates the likelihoods. */
BiobaseLikelihoods
create_likelihoods_from_counts(const BiobaseCounts & counts);


/** Given a vector of likelihoods, calculates the cumulative likelihoods... I.e. the likelihoods
of a particular score or better. */
BiobaseLikelihoods
create_cumulative_likelihoods_from_non(const BiobaseLikelihoods & likelihoods);


/** Given a vector of cumulative counts, calculates the likelihoods.. */
BiobaseLikelihoods
create_likelihoods_from_cumulative_counts(const BiobaseCounts & counts);


float_t
get_likelihood(const BiobaseLikelihoods & likelihoods, float_t score);


/** Functor to get probability of particular score from likelihoods. */
struct QuantisedScores
{
	const BiobaseLikelihoods * likelihoods;

	QuantisedScores( const BiobaseLikelihoods * likelihoods = 0 );

	float_t operator()( float_t score ) const;
};


/** Get quantised score functor for particular pssm. */
QuantisedScores
get_biobase_quantised_scores(
	const TableLink & key,
	bool background,
	bool or_better );



/** Stores and persists quantised counts for pssms indexed by key. Can generate Ott normalisations and likelihoods
generated from these counts on demand. */
class LikelihoodsCache
	: public Singleton< LikelihoodsCache >
{
	friend struct Singleton< LikelihoodsCache >;

public:
	typedef boost::shared_ptr< LikelihoodsCache > ptr_t;
	typedef TableLink key_t;

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & counts;
		//ott_normalisers & likelihoods are generated on demand and do not need to be persisted
    }

	void init_singleton();

protected:
	typedef std::map< key_t, BiobaseCounts > count_map_t;
	typedef std::map< key_t, BiobaseLikelihoods > likelihood_map_t;

	count_map_t counts;
	likelihood_map_t background_likelihoods;
	likelihood_map_t background_likelihoods_or_better;
	likelihood_map_t binding_likelihoods;
	likelihood_map_t binding_likelihoods_or_better;

public:

	/** Quantise counts for all pssms in the iterators. */
	template <class PssmIt>
	void
	quantise_counts(
		PssmIt pssm_begin,
		PssmIt pssm_end,
		const seq_t & seq);

	/** Gets the counts for the given key.
	The counts are the number of times particular scores are seen when testing the pssm against a
	random background sequence.
	Creates an empty vector of 0 counts if not found. 
	Calling this invalidates any stored ott normalisations or likelihoods.
	These will be recalculated as needed. */
	BiobaseCounts *
	get_counts( const key_t & key );

	/**
	Gets the total counts for the given key.
	*/
	unsigned
	get_total_counts( const key_t & key ) const;

	/** Gets the likelihoods of particular scores in a background random or a binding sequence for the named pssm.
	 * Calculates them if needed. Returns 0 if cannot.
	 * @param key Name of the pssm */
	const BiobaseLikelihoods *
	get_score_likelihoods(
		const key_t & key,
		bool background,
		bool or_better );

	/** Gets the likelihoods of particular scores or better in a background random sequence for the named pssm.
	Calculates them from counts if needed. Returns 0 if cannot. */
	const BiobaseLikelihoods *
	get_background_score_or_better_likelihoods( const key_t & key );

	/** Gets the likelihoods of biobase scores assuming a sequence matches the consensus of a named pssm. */
	const BiobaseLikelihoods *
	get_binding_score_likelihoods( const key_t & key );

	/** Gets the likelihoods of biobase scores (or better) assuming a sequence matches the consensus of a named pssm. */
	const BiobaseLikelihoods *
	get_binding_score_or_better_likelihoods( const key_t & key );

	/** Updates the counts of pssms in the biobase db. */
	void update_counts( BiobaseDb & db, const seq_t & seq );

	bool operator==( const LikelihoodsCache & rhs ) const;
};




/** Quantises and counts the scores a pssm has over a sequence. Returns the number of scores generated. */
template <class PssmT, class SeqIt>
size_t
quantise_scores(
	const PssmT & pssm,
	SeqIt seq_begin,
	SeqIt seq_end,
	BiobaseCounts & quantised_counts)
{
	const size_t num_quanta = quantised_counts.size();
	if (0 == num_quanta)
	{
		throw std::logic_error( "0 == num_quanta" );
	}

	size_t num_samples = 0;

	//score each sub-sequence and its complement to count frequencies of each quantum
	for (; (size_t) (seq_end - seq_begin) >= pssm.size(); ++seq_begin) {

		//match uncomplemented
		{
			const float_t score = pssm.score(seq_begin, false);
			const size_t idx = get_biobase_score_index(num_quanta, score);
			quantised_counts[idx]++; //increment the count for the indexed bucket
			++num_samples;
		}

		//match complemented
		{
			const float_t score = pssm.score(seq_begin, true);
			const size_t idx = get_biobase_score_index(num_quanta, score);
			quantised_counts[idx]++; //increment the count for the indexed bucket
			++num_samples;
		}
	}

	return num_samples;
}



template <class PssmIt>
void
LikelihoodsCache::quantise_counts(
	PssmIt pssm_begin,
	PssmIt pssm_end,
	const seq_t & seq)
{
	//for each pssm
	for ( ; pssm_begin != pssm_end; ++pssm_begin)
	{
		BiobaseCounts * counts = get_counts( pssm_begin->first );

		quantise_scores( make_pssm(pssm_begin->second), seq.begin(), seq.end(), *counts );
	}
}



BIO_NS_END

#endif //BIO_BIOBASE_LIKELIHOOD_H_
