/**
@file

Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"

#include "bio/application.h"
#include "bio/biobase_filter.h"
USING_BIO_NS;

#include "biopsy/analyse.h"
#include "biopsy/transfac.h"
using namespace biopsy;


using namespace boost;
using namespace std;

double add_log_likelihoods( double ll1, double ll2 )
{
	return log( exp( ll1 ) + exp( ll2 ) );
}


typedef std::map< double, double > likelihood_map; /**< Maps scores to likelihoods. */
typedef likelihood_map::value_type likelihood_entry;
typedef boost::shared_ptr< likelihood_map > likelihood_map_ptr;

void insert_score_into_likelihood_map( double score, double log_likelihood, likelihood_map & c )
{
	typedef likelihood_map::iterator iterator;
	typedef likelihood_map::value_type entry;
	iterator i = c.lower_bound( score );
	if( c.end() == i || score < i->first ) {
		//we have a new likelihood to insert
		c.insert( entry( score, log_likelihood ) );
	} else {
		//update the existing one (which has the same score)
		i->second = add_log_likelihoods( log_likelihood, i->second );
	}
}

void reduce_likelihood_map_size( likelihood_map & c, unsigned max_size )
{
	int num_to_remove = c.size() - max_size;
	while( num_to_remove > 0 )
	{
		//what range do we have?
		const double max_score = c.rbegin()->first;
		const double min_score = c.begin()->first;
		BOOST_ASSERT( min_score <= max_score );

		//we will consolidate all consecutive scores whose difference is smaller than this quantum
		const double quantum = (max_score - min_score) / max_size;

		//find the differences between consecutive iterators
		likelihood_map::iterator i = c.begin();
		while( num_to_remove > 0 )
		{
			likelihood_map::iterator j = i;
			++j;
			if( c.end() == j ) break;

			const double diff = j->first - i->first;
			BOOST_ASSERT( diff > 0. );
			if( diff >= quantum )
			{
				//nothing to do...
				++i;
			}
			else
			{
				//had to be careful here with floating point underflow
				const double score1 = i->first;
				const double score2 = j->first;
				BOOST_ASSERT(score1 != score2);
				const double prob1 = exp(i->second);
				const double prob2 = exp(j->second);
				const double prob_sum = prob1 + prob2;
				const double relative_prob2 = prob2 / prob_sum;
				const double weighted_avg_score = score1 + (score2 - score1) * relative_prob2;
				BOOST_ASSERT(score1 <= weighted_avg_score);
				BOOST_ASSERT(weighted_avg_score <= score2);

				c.erase( i );
				c.erase( j );

				std::pair< likelihood_map::iterator, bool > insert_result = c.insert( likelihood_entry( weighted_avg_score, log( prob_sum ) ) );
				BOOST_ASSERT( insert_result.second );

				--num_to_remove;
				i = insert_result.first;
			}
		}
	}
}

template< typename OutputIt >
void
get_pssm_likelihoods(
    pssm_ptr _pssm,
    unsigned num_output_bins,
    OutputIt output_it,
    unsigned max_likelihood_map_size = 50000
) {
	//initially we have likelihood 1 of seeing score 0.
	likelihood_map_ptr last_scores( new likelihood_map );
	last_scores->insert( likelihood_entry( 0., 0. ) );

	likelihood_map_ptr new_scores( new likelihood_map );

	const double log_quarter = log( .25 );

	BOOST_FOREACH( const nucleo_dist & dist, *_pssm )
	{
		//boost::timer _timer;
		new_scores->clear();

		BOOST_FOREACH( likelihood_entry last, *last_scores )
		{
			for( unsigned i = 0; 4 != i; ++i )
			{
				insert_score_into_likelihood_map( last.first + dist.get(i), last.second + log_quarter, *new_scores );
			}
		}

		//cout << "Have got " << new_scores->size() << " scores, took " << _timer.elapsed() << " secs.\n";
		//_timer.restart();
		reduce_likelihood_map_size( *new_scores, max_likelihood_map_size );
		//cout << "Reduced to " << new_scores->size() << " scores, took " << _timer.elapsed() << " secs.\n";
		new_scores.swap( last_scores );
	}
}



/**
Tests how long it takes to calculate likelihood of PSSM scores under a background model.
*/
struct TestPssmLikelihoodsApp : Application
{
	int task()
	{
		string_vec_ptr pssm_names = get_transfac_pssm_accessions(BiobasePssmFilter::get_all_pssms_filter());
		BOOST_FOREACH( string name, *pssm_names )
		{
			const pssm_info & pssm_info = get_pssm( name );
			cout << name;

			std::vector< double > likelihoods;
			boost::timer _timer;
			get_pssm_likelihoods( pssm_info._pssm, 100, back_inserter( likelihoods ), 5000 );
			cout << " took " << _timer.elapsed() << " secs.\n";
		}

		return 0;
	}
};





int
main(
	int argc,
	char * argv[] )
{
	return TestPssmLikelihoodsApp().main( argc, argv );
}



