/* Copyright John Reid 2007, 2011
*/

#include "bio-pch.h"


#include "bio/defs.h"
//#define BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP


#include "bio/biobase_likelihoods.h"
#include "bio/environment.h"
#include "bio/biobase_match.h"
#include "bio/biobase_db.h"
#include "bio/biobase_filter.h"
#include "bio/serialisable.h"
#include "bio/pssm_likelihood_cache.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/filesystem/path.hpp>
using namespace boost;
namespace fs = boost::filesystem;

#include <fstream>
#include <numeric>
#include <iostream>
using namespace std;

#if BOOST_VERSION >= 104400
#define _FPC_NS ::boost::math::fpc
#else
# define _FPC_NS ::boost::test_tools
#endif


BIO_NS_START

/** Gets the index of a given score in [0,1]. */
size_t get_biobase_score_index(size_t size, float_t score)
{
	if (size == 0)
	{
		throw std::logic_error( "size == 0" );
	}

	if (0.0 > score)
	{
		throw std::logic_error( BIO_MAKE_STRING( "Score < 0: " << score ) );
	}

	if (score * 2 * size > 2 * size + 1)
	{
		throw std::logic_error( BIO_MAKE_STRING( "Score > 1 and a bit: " << score ) );
	}

	size_t idx = (size_t) (score * size);
	if (size == idx) //cater for a perfect match
	{
		--idx;
	}
	assert(0 <= idx && idx < size);

	return idx;
}

float_t
get_likelihood(const BiobaseLikelihoods & likelihoods, float_t score)
{
	return likelihoods[get_biobase_score_index(likelihoods.size(), score)];
}



/** Turn the counts of quantised scores into cumulative counts. */
BiobaseCounts
create_cumulative_from_counts(const BiobaseCounts & counts)
{
	BiobaseCounts result;
	size_t cumulative = 0;
	for (BiobaseCounts::const_iterator i = counts.begin();
		i != counts.end();
		++i)
	{
		cumulative += *i;
		result.push_back(cumulative);
	}
	return result;
}


/** Given a vector of counts, calculates the likelihoods. */
BiobaseLikelihoods
create_likelihoods_from_counts(const BiobaseCounts & counts)
{
	const size_t num_samples = std::accumulate(counts.begin(), counts.end(), 0);
	if (0 == num_samples)
	{
		throw std::logic_error( "No samples to generate likelihoods from" );
	}

	BiobaseLikelihoods result;
	for (BiobaseCounts::const_iterator i = counts.begin();
		counts.end() != i;
		++i)
	{
		result.push_back(((float_t) *i) / ((float_t) num_samples));
	}

	return result;
}


BiobaseLikelihoods
create_cumulative_likelihoods_from_non(const BiobaseLikelihoods & likelihoods)
{

	float_t cumulative = 0.0;
	unsigned idx = likelihoods.size();
	BiobaseLikelihoods result( idx );
	for (BiobaseLikelihoods::const_reverse_iterator i = likelihoods.rbegin();
		likelihoods.rend() != i;
		++i)
	{
		cumulative += *i;
		//adjust for rounding error
		cumulative = std::min( float_t( 1.0 ), cumulative );

		result[ --idx ] = cumulative;
	}
	BOOST_ASSERT(
		boost::test_tools::check_is_close(
			float_t( 1.0 ),
			cumulative,
			BIO_FPC_NS::percent_tolerance( 0.01f ) ) );

	return result;
}



unsigned
LikelihoodsCache::get_total_counts( const key_t & key ) const
{
	//look for the counts
	count_map_t::const_iterator c = counts.find( key );
	if ( counts.end() == c )
	{
		//we could not find the counts
		return 0;
	}

	//counts already in map
	return std::accumulate( c->second.begin(), c->second.end(), 0 );
}



BiobaseCounts *
LikelihoodsCache::get_counts(const key_t & key)
{
	BiobaseCounts * result = 0;

	//look for the counts
	count_map_t::iterator c = counts.find(key);
	if (counts.end() == c)
	{
		//we could not find the counts - insert a vector of 0's
		result =
			&(counts.insert(
				make_pair(
					key,
					BiobaseCounts(BioEnvironment::singleton().num_normalisation_quanta, 0))).first->second);
	}
	else
	{
		//counts already in map
		result = &(c->second); //retrieve them from map
	}

	//invalidate likelihoods
	background_likelihoods.erase(key);
	background_likelihoods_or_better.erase(key);
	binding_likelihoods.erase(key);
	binding_likelihoods_or_better.erase(key);

	return result;
}




const BiobaseLikelihoods *
LikelihoodsCache::get_score_likelihoods(
	const key_t & key,
	bool background,
	bool or_better)
{
	//which map are we going to look for the likelihoods in
	likelihood_map_t & map =
		background
			? (or_better ? background_likelihoods_or_better : background_likelihoods)
			: (or_better ? binding_likelihoods_or_better : binding_likelihoods);

	const BiobaseLikelihoods * result = 0;

	//look for the likelihoods
	likelihood_map_t::const_iterator l = map.find(key);
	if (map.end() == l)
	{
		//we could not find the likelihoods
		BiobaseLikelihoods likelihoods;

		if (or_better)
		{
			//we need cumulative likelihoods

			//first get the non-cumulative
			const BiobaseLikelihoods * non_cumulative = get_score_likelihoods(key, background, false);
			if (0 != non_cumulative)
			{
				//calculate the cumulative from the non.
				likelihoods = create_cumulative_likelihoods_from_non(*non_cumulative);
			}
		}
		else if (background)
		{
			//can we generate the likelihoods?
			count_map_t::const_iterator c = counts.find(key);
			if (counts.end() != c)
			{
				likelihoods = create_likelihoods_from_counts(c->second);
			}
		}
		else
		{
			likelihoods = *PssmLikelihoodCache::singleton().get_likelihoods(key);
		}

		//did we calculate any likelihoods?
		if (likelihoods.size() > 0)
		{
			//insert them
			std::pair<likelihood_map_t::iterator, bool> insert_result = map.insert(std::make_pair(key, likelihoods));
			assert(insert_result.second); //make sure it wasn't already in the map
			result = &(insert_result.first->second);
		}
	}
	else
	{
		//likelihoods already in map
		result = &(l->second); //retrieve normalisation from map
	}

	return result;
}

void
LikelihoodsCache::init_singleton()
{
	try
	{
		deserialise< false >(
			*this,
			fs::path(
				BioEnvironment::singleton().get_likelihoods_cache_file()
			)
		);
	}
	catch( const std::exception & ex )
	{
		std::cout << "LikelihoodsCache::init_singleton(): could not deserialise: " << ex.what() << std::endl;
	}
	catch( ... )
	{
		std::cout << "LikelihoodsCache::init_singleton(): could not deserialise: unknown error" << std::endl;
	}
}


void
LikelihoodsCache::update_counts(BiobaseDb & db, const seq_t & seq)
{
	//Normalising matrices
	quantise_counts(
		matrix_filter_it(db.get_matrices().begin(), db.get_matrices().end()),
		matrix_filter_it(db.get_matrices().end(), db.get_matrices().end()),
		seq);

	//Normalising sites
	quantise_counts(
		site_filter_it(db.get_sites().begin(), db.get_sites().end()),
		site_filter_it(db.get_sites().end(), db.get_sites().end()),
		seq);
}

bool
LikelihoodsCache::operator==(const LikelihoodsCache & rhs) const
{
	return counts == rhs.counts;
}




QuantisedScores::QuantisedScores( const BiobaseLikelihoods * likelihoods )
	: likelihoods( likelihoods )
{
}


float_t
QuantisedScores::operator()( float_t score ) const
{
	if( 0 == likelihoods )
	{
		throw std::logic_error( "Null pointer in QuantisedScores::operator()" );
	}

	return get_likelihood( *likelihoods, score );
}

QuantisedScores
get_biobase_quantised_scores(
	const LikelihoodsCache::key_t & key,
	bool background,
	bool or_better )
{
	return QuantisedScores( LikelihoodsCache::singleton().get_score_likelihoods( key, background, or_better ) );
}



BIO_NS_END

