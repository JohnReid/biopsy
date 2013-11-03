#ifndef BIO_PAIR_ANALYSIS_H_
#define BIO_PAIR_ANALYSIS_H_

#include "bio/defs.h"
#include "bio/counter.h"
#include "bio/transcription_factor.h"
#include "bio/score_map.h"
#include "bio/iterator.h"
#include "bio/math.h"

#include <boost/operators.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>


BIO_NS_START





/** Contains statistics about the distribution of pairs of binders. */
template< typename Binder >
struct PairStatistics
{
	typedef Binder binder_t;
	typedef typename BindingHitSet< binder_t >::type hit_set_t;
	typedef boost::tuple< const binder_t *, const binder_t * > binder_pair_t;
	typedef BindingHit< binder_t > hit_t;
	typedef Counter< const binder_t * > binder_counter_t; /**< Counts instances of a single factor. */
	typedef Counter< binder_pair_t > pair_counter_t; /**< Counts instances of pairs of factors. */
	typedef Counter< unsigned > distance_counter_t; /**< Counts instances of distances between a pair of factors. */
	typedef std::map< binder_pair_t, distance_counter_t > distances_map_t; /**< Stores many distance counts. */
	typedef Counter< std::string > sequence_counter_t; /**< A counter that counts how many hits in each sequence (indexed by name of sequence). */
	typedef std::map< const binder_t * , sequence_counter_t > sequence_counts_map_t; /**< Keep a count for each pssm. */
	typedef std::map< std::string, unsigned > seq_lengths_map_t; /**< Remembers the length of sequences. */
	typedef typename ScoreMap< binder_pair_t >::type pair_evidences_map_t; /**< Used to rank findings in order of interest. */
	typedef typename ScoreMap< const binder_t * >::type single_evidences_map_t; /**< Used to rank findings in order of interest. */

	binder_counter_t binder_counts;
	pair_counter_t pair_counts;
	distances_map_t distances_map;
	sequence_counts_map_t sequence_counts;
	seq_lengths_map_t seq_lengths;
	unsigned distance;
	unsigned num_remos;
	unsigned num_bases;

	PairStatistics( unsigned distance = 0, unsigned num_remos = 0, unsigned num_bases = 0 );

	/** The likelihood that 2 hits are within the distance assuming a uniform distribution. */
	double get_p() const;

	/** Add the hits to the counts. */
	template< 
		typename HitIt
	>
	void add_hits(
		const std::string & sequence_group_name,
		HitIt results_begin, //must be ordered by position
		HitIt results_end,
		unsigned seq_length );

	double get_pair_log_likelihood_ratio( unsigned pair_count, unsigned f1_count, unsigned f2_count, double p );

	void analyse( unsigned num_output );
	void analyse_distances( unsigned num_output );
	void analyse_clusters( unsigned num_output );
	void analyse_singles( unsigned num_output );

protected:
	friend class boost::serialization::access;
    template< typename Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
		throw std::logic_error( "PairStatistics::serialize() not implemented" );
		///**
        //ar & binder_counts;
		//ar & pair_counts;
		//ar & distances_map;
		//ar & sequence_counts;
		ar & seq_lengths;
		ar & distance;
		ar & num_remos;
		ar & num_bases;
		//*/
    }
};





template< typename Binder >
PairStatistics< Binder >::PairStatistics(unsigned distance, unsigned num_remos, unsigned num_bases)
	: distance(distance)
	, num_remos(num_remos)
	, num_bases(num_bases)
{
}

/** The likelihood that 2 hits are within the distance assuming a uniform distribution. */
template< typename Binder >
double PairStatistics< Binder >::get_p() const
{
	const double numerator = distance;
	const double denominator = (double(num_bases) - double(num_remos) * double(distance) / 2.0);

	if (0.0 == denominator)
	{
		throw std::logic_error("Divide by 0");
	}

	const double result = numerator / denominator;
	if (! BIO_FINITE(result))
	{
		std::cout
			<< "distance=" << distance
			<< ", num_remos=" << num_remos
			<< ", num_bases=" << num_bases
			<< "\n"; 
		throw std::logic_error("Divide by 0");
	}

	return result;
}

/** Add the hits to the counts. */
template< 
	typename Binder
>
template< 
	typename HitIt
>
void PairStatistics< Binder >::add_hits(
	const std::string & sequence_group_name,
	HitIt results_begin, //must be ordered by position
	HitIt results_end,
	unsigned seq_length)
{
	//add or initialise the sequence length
	if (seq_lengths.end() == seq_lengths.find(sequence_group_name))
	{
		seq_lengths[sequence_group_name] = seq_length;
	}
	else
	{
		seq_lengths[sequence_group_name] += seq_length;
	}

	//for each hit
	for( HitIt h1 = results_begin; results_end != h1; ++h1 )
	{
		//update the sequence counts
		sequence_counts[ h1->get_binder() ].increment( sequence_group_name );

		binder_counts.increment( h1->get_binder() );

		//look for other hits ahead of h1 but inside the distance
		for( HitIt h2 = boost::next( h1 );
			results_end != h2;
			++h2)
		{
			//check h2 is ahead of end of h1
			if( h2->get_position() < h1->get_end() )
			{
				//it isn't
				continue;
			}

			//check h2 is not too far ahead, i.e. inside the distance
			if( h2->get_position() >= h1->get_end() + int( distance ) )
			{
				//it is too far ahead
				break;
			}

			//count them
			const binder_pair_t pair( h1->get_binder(), h2->get_binder() );
			pair_counts.increment( pair );

			//are we analysing the distances and we are not too close to the start or end of the sequence
			//to skew the distribution?
			if ( h2->get_position() >= int( distance ) && h1->get_end() + distance < seq_length)
			{
				distances_map[ pair ].increment( h2->get_position() - ( h1->get_position() + h1->get_length() ) );
			}
		}
	}
}

/** Get log likelihood ratio for counts of pair of transcription factors. */
template< typename Binder >
double PairStatistics< Binder >::get_pair_log_likelihood_ratio(unsigned pair_count, unsigned f1_count, unsigned f2_count, double p)
{
	//cout << pair_count << "\t" << f1_count << "\t" << f2_count << "\t" << p << "\n";

	const std::vector< unsigned > obs = boost::assign::list_of( pair_count )( f1_count * f2_count );
	const std::vector< double > ps = boost::assign::list_of( p )( 1.0 - p );

	return
		calc_ln_multinomial_evidence(
			obs.begin(),
			obs.end(),
			ps.begin(),
			ps.begin() );

}

template< typename Binder >
void PairStatistics< Binder >::analyse(unsigned num_output)
{
	analyse_distances(num_output);
	analyse_singles(num_output);
	analyse_clusters(num_output);
}

template< typename Binder >
void PairStatistics< Binder >::analyse_singles(unsigned num_output)
{
	if (0 == num_bases)
	{
		throw std::logic_error("num_bases not set");
	}

	//generate one alpha for all pssms
	std::vector< double > alpha;
	for (seq_lengths_map_t::const_iterator l = seq_lengths.begin();
		seq_lengths.end() != l;
		++l)
	{
		alpha.push_back(double(l->second) / double(num_bases));
	}

	//std::cout << "Calculating scores...\n";
	single_evidences_map_t scores;
	for (typename sequence_counts_map_t::const_iterator p = sequence_counts.begin();
		sequence_counts.end() != p;
		++p)
	{
		const sequence_counter_t & counts = p->second;

		//copy the counts into the observation vector
		std::vector< unsigned > obs;
		for (seq_lengths_map_t::const_iterator s = seq_lengths.begin();
			seq_lengths.end() != s;
			++s)
		{
			sequence_counter_t::const_iterator c = counts.find(s->first);
			obs.push_back(counts.end() == c ? 0 : c->second);
		}

		const double uniform_evidence = calc_ln_multinomial_evidence(obs.begin(), obs.end(), alpha.begin(), alpha.begin());
		scores.insert( single_evidences_map_t::value_type( p->first, uniform_evidence ));
	}

	std::cout << "Printing sequence scores...\n";
	unsigned i = 0;
	BOOST_FOREACH(
		const typename ScoreMap< const binder_t * >::value_type & evidence,
		scores.get< ScoreMap< const binder_t * >::score >() )
	{
		std::cout << evidence.score << "\t" << evidence.key << "\n";
		sequence_counts[ evidence.key ].print( true, std::cout, "Sequence", false, 20, 0, true );
		std::cout << "\n";

		++i;

		//only print so many
		if( num_output == i )
		{
			break;
		}
	}
}

template< typename Binder >
void PairStatistics< Binder >::analyse_clusters(unsigned num_output)
{
	const double p = get_p();
	std::cout << "p = " << p << "\n";

	//look for evidence that they cluster together...
	pair_evidences_map_t evidences;
	for (pair_counter_t::const_iterator c = pair_counts.begin();
		pair_counts.end() != c;
		++c)
	{
		const unsigned count_1 = binder_counts.get_count( c->first.get< 0 >() );
		const unsigned count_2 = binder_counts.get_count( c->first.get< 1 >() );
		const unsigned pair_count = c->second;

		const double evidence =
			get_pair_log_likelihood_ratio(
				pair_count,
				count_1,
				count_2,
				p);

		evidences.insert( pair_evidences_map_t::value_type( c->first, evidence ) );
	}

	//print the most interesting...
	unsigned i = 0;
	for (pair_evidences_map_t::const_iterator e = evidences.begin();
		evidences.end() != e && num_output != i;
		++e, ++i)
	{
		std::cout
			<< e->score
			<< "\t" << pair_counts.get_count( e->key )
			<< "\t" << binder_counts.get_count( e->key.get< 0 >() )
			<< "\t" << binder_counts.get_count( e->key.get< 1 >() )
			<< "\t" << e->key
			<< "\n";
	}
}

template< typename Binder >
void PairStatistics< Binder >::analyse_distances( unsigned num_output )
{
	pair_evidences_map_t evidences;
	double alpha = 1.0;

	for (distances_map_t::const_iterator d = distances_map.begin();
		distances_map.end() != d;
		++d)
	{
		const distance_counter_t & distance_counter = d->second;

		std::vector< unsigned > distance_counts;
		for (unsigned i = 0; distance != i; ++i)
		{
			distance_counts.push_back(distance_counter.get_count(i));
		}
		const double uniform_evidence =
			calc_ln_multinomial_uniform_evidence(
				distance_counts.begin(),
				distance_counts.end(),
				single_value_iterator< double >(alpha));

		evidences.insert( pair_evidences_map_t::value_type( d->first, uniform_evidence ) );
	}

	//print the most interesting...
	unsigned i = 0;
	BOOST_FOREACH( const ScoreMap< binder_pair_t >::value_type & evidence, evidences.get< ScoreMap< binder_pair_t >::score >() )
	{
		std::cout << distances_map[ evidence.key ].get_total() << " " << evidence.score << " " << evidence.key << "\n";
		for( unsigned j = 0; distance != j; ++j )
		{
			distances_map[ evidence.key ].initialise( j );
		}
		distances_map[ evidence.key ].print( false, std::cout, "Distance", false, 0, 60 );
		std::cout << "\n";

		++i;

		//only print so many
		if( num_output == i )
		{
			break;
		}
	}
}



BIO_NS_END

#endif //BIO_PAIR_ANALYSIS_H_

