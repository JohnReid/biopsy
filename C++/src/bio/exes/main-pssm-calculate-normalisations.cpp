/**
@file

Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"



#include <bio/application.h>
#include <bio/hmm_gen_sequence.h>
#include <bio/hmm_dna.h>
#include <bio/options.h>
#include <bio/random.h>
#include <bio/biobase_filter.h>
#include <bio/environment.h>
#include <bio/biobase_likelihoods.h>
#include <bio/serialisable.h>
USING_BIO_NS;


#include <boost/progress.hpp>
#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


//#include <windows.h> 

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

/**
 * Calculates normalisations that are used to score PSSMs using old(ish) scoring method. Generates random sequences and then
 * applies PSSMs to them in order to empirically estimate a distribution of scores (or likeilhoods?).
 * Can either generate sequences from a uniform background distribution or from higher-order Markov models trained
 * on specific species. Will run indefinitely until interrupted with Ctrl-C and will store estimates on a regular basis or when
 * requested to using Ctrl-BREAK.
 */
struct CalculateNormalisationsApp : Application
{
	size_t seq_length;
	bool use_uniform_dist;
	unsigned hmm_num_states;
	unsigned hmm_order;
	unsigned serialise_every_so_often;
	bool want_to_exit;
	bool want_to_serialise;
	bool have_reported_count_mismatch;

	CalculateNormalisationsApp()
		: want_to_exit( false )
		, want_to_serialise( false )
		, have_reported_count_mismatch( false )
	{
		get_options().add_options()
			("seq_length", po::value(&seq_length)->default_value(2000), "length of sequence to normalise over")
			("use_uniform_dist", po::bool_switch(&use_uniform_dist)->default_value(false), "use a uniform distribution")
			("hmm_num_states", po::value(&hmm_num_states)->default_value(1), "number of states in HMM")
			("hmm_order", po::value(&hmm_order)->default_value(3), "order of HMM")
			("serialise,s", po::value(&serialise_every_so_often)->default_value(0), "serialise every so often (s)")
			;
	}

	virtual bool ctrl_handler( CtrlSignal signal )
	{
		want_to_serialise = true;
		want_to_exit = CTRL_BREAK_SIGNAL != signal;

		if ( want_to_exit )
		{
			cout << "Will store values and exit at next iteration" << endl;
		}
		else
		{
			cout << "Will store values and continue at next iteration" << endl;
		}

		return true;
	}

	struct has_total_counts_at_most
	{
		unsigned num;
		has_total_counts_at_most( unsigned num ) : num( num ) { }
		template< typename MapValue >
		bool operator()( const MapValue & value ) const {
			return num >= LikelihoodsCache::singleton().get_total_counts( value.first );
		}
	};

	template< typename PssmIt >
	void calculate_scores( 
		PssmIt pssms_begin,
		PssmIt pssms_end,
		const seq_t & norm_seq)
	{
		//first see what the maximum and minimum # scores for each pssm
		unsigned min_total_counts = pssms_end == pssms_begin ? 0 : std::numeric_limits< unsigned >::max();
		unsigned max_total_counts = 0;
		for( PssmIt p = pssms_begin;
			pssms_end != p;
			++p )
		{
			const TableLink link = p->first;
			const unsigned total_counts = LikelihoodsCache::singleton().get_total_counts( link );

			min_total_counts = std::min( total_counts, min_total_counts );
			max_total_counts = std::max( total_counts, max_total_counts );
		}

		if( min_total_counts != max_total_counts && ! have_reported_count_mismatch )
		{
			std::cout
				<< "Smallest # normalisation counts: " << min_total_counts << endl
				<< "Largest # normalisation counts: " << max_total_counts << endl;

			have_reported_count_mismatch = true;
		}

		//restrict the iterators to those with the counts at most a weighted avg of min and max
		const unsigned count_threshold = ( min_total_counts + 9 * max_total_counts ) / 10;

		LikelihoodsCache::singleton().quantise_counts(
			make_filter_iterator( has_total_counts_at_most( count_threshold ), pssms_begin, pssms_end ),
			make_filter_iterator( has_total_counts_at_most( count_threshold ), pssms_end, pssms_end ),
			norm_seq );
	}

	int task()
	{
		cout << "Using sequence length of " << seq_length << endl;
		if (use_uniform_dist)
		{
			cout << "Will generate sequences from uniform random distribution" << endl;
		}
		else
		{
			cout << "Will generate sequences from species HMMs" << endl;
		}


		if( serialise_every_so_often > 0 )
		{
			cout << "Will serialise normalisations every " << serialise_every_so_often << " seconds\n";
		}
		else
		{
			cout << "Will only serialise normalisations on request\n";
		}

		register_ctrl_handler();
		cout
			<< endl
			<< "Hit Ctrl-BREAK to save current state" << endl
			<< "Hit Ctrl-C to save current state and exit" << endl
			<< endl;

		const BiobasePssmFilter filter = BiobasePssmFilter::get_all_pssms_filter();

		//use the timer to decide whether to serialise every so often
		boost::timer timer;

		//repeat until user breaks
		cout << "Updating counts over random sequences" << endl;
		while( true )
		{
			seq_t norm_seq;
			norm_seq.reserve( seq_length );
			{
				if ( use_uniform_dist )
				{
					generate_random_nucleotide_seq( inserter( norm_seq, norm_seq.begin() ), seq_length );
				}
				else
				{
					DnaHmmOrderNumStateMap & hmm_map = DnaHmmOrderNumStateMap::singleton();
					if (! hmm_map.contains_model(hmm_num_states, hmm_order))
					{
						throw std::logic_error( "Do not have model with that # states and order" );
					}
					hmm_map.get_model( hmm_num_states, hmm_order ).append_random_sequence( norm_seq, seq_length );
				}
			}

			//check we have some pssms to work on
			if( get_matrices_begin( filter ) == get_matrices_end( filter )
				&& get_sites_begin( filter ) == get_sites_end( filter ) )
			{
				throw std::logic_error( "No pssms to work on - pssm filter too restrictive" );
			}

			//estimate the scores for each matrix and site on the sequence
			calculate_scores(
				get_matrices_begin( filter ),
				get_matrices_end( filter ),
				norm_seq);
			calculate_scores(
				get_sites_begin( filter ),
				get_sites_end( filter ),
				norm_seq);

			//do we want to serialise because we have been running for so long?
			if( serialise_every_so_often > 0 && timer.elapsed() > double( serialise_every_so_often ) )
			{
				want_to_serialise = true;

				timer.restart();
			}

			//are we going to serialise the scores?
			if( want_to_serialise )
			{
				cout << "Storing values" << endl;
				serialise< false >(
					LikelihoodsCache::singleton(),
					boost::filesystem::path(
						BioEnvironment::singleton().get_likelihoods_cache_file()
					)
				);

				want_to_serialise = false;
				have_reported_count_mismatch = false;

				//do we want to exit?
				if( want_to_exit )
				{
					break;
				}
			}
		}

		return 0;
	}
};


int
main(int argc, char * argv [])
{
	return CalculateNormalisationsApp().main(argc, argv);
}
