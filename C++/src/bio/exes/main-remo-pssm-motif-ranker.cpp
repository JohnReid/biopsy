/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/pssm_motif.h"
#include "bio/remo_analysis.h"
USING_BIO_NS

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;



struct PssmMotifRankerApp : Application, AnalysisVisitor
{
	typedef std::vector< std::string > string_vec_t;
	typedef std::multimap< BIO_NS::float_t, std::string > rank_map_t; /** Maps scores to sequence or remo names. */
	typedef boost::tuples::tuple< rank_map_t, rank_map_t > motif_ranking_t; /** Contains scores for sequences and remos. Sequences first. */
	typedef std::map< PssmMotif::set_t, motif_ranking_t > rank_map_collection_t; /** Maps a motif set to its ranks. */
	typedef std::set< PssmMotif::set_t > motif_set_set_t; /** A collection of pssm motif sets we're interested in. */
	typedef std::map< PssmMotif::ptr_t, std::string > motif_desc_map_t; /** Maps motifs to their descriptions. */

	//members defined by program arguments
	string_vec_t motif_descriptions;
	unsigned num_to_display;
	double prior;
	bool leave_one_out;

	PssmMotif::vec_t pssm_motifs;
	rank_map_collection_t rank_map;
	PssmMotif::score_map_t sequence_evidence;
	motif_set_set_t motif_sets;
	motif_desc_map_t motif_descriptions_map;

	PssmMotifRankerApp()
	{
		get_options().add_options()
			("motif,m", po::value(&motif_descriptions), "motif description")
			("num_to_display,n", po::value(&num_to_display)->default_value(10), "# to display")
			("prior,p", po::value(&prior)->default_value(100.0), "prior")
			("leave_one_out,l", po::value(&leave_one_out)->default_value(true), "provides rankings for subsets of the motifs")
			//("open_ensembl,e", po::value(&open_ensembl)->default_value(false), "open ensembl webpages in firefox")
			;

		add_analysis_options(get_options());

		get_positional_options().add("motif", -1);
	}



	void parse_motif_descriptions()
	{
		for (string_vec_t::const_iterator m = motif_descriptions.begin();
			motif_descriptions.end() != m;
			++m)
		{
			PssmMotif::ptr_t motif = PssmMotif::parse(*m);
			pssm_motifs.push_back(motif);
			motif_descriptions_map[motif] = *m;
		}
	}



	bool visit_sequence_group(const std::string & seq_group_name)
	{
		sequence_evidence.clear();

		return true;
	}



	void leave_sequence_group(const std::string & seq_group_name)
	{
		add_score_to_ranks(sequence_evidence, true, seq_group_name);
	}




	void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		match_result_vec_t & results,
		const seq_t & sequence)
	{
		PssmMotif::score_map_t remo_evidence;

		sort_by_position(results);

		typedef std::set<MatchResults> hit_set_t;
		hit_set_t hits_in_motifs_set;

		//for each motif
		BIO_NS::float_t min_score = 1.0f;
		for (PssmMotif::vec_t::const_iterator m = pssm_motifs.begin();
			pssm_motifs.end() != m;
			++m)
		{
			Score & remo_score = remo_evidence[*m];
			Score & sequence_score = sequence_evidence[*m];

			//look for the motifs
			PssmMotif::HitVec hits;
			m->get()->find_in(results, hits);

			//for each hit
			for (PssmMotif::HitVec::const_iterator hit = hits.begin();
				hits.end() != hit;
				++hit)
			{
				const double hit_score = double(PssmMotif::score(*hit));
				const double evidence = hit_score / prior;

				sequence_score.add(evidence);
				remo_score.add(evidence);
			}
		}

		add_score_to_ranks(remo_evidence, false, remo_name);

	}



	/** Create a collection of motif sets we're interested in. */
	void build_interesting_motif_sets()
	{
		motif_sets.clear();

		//We're always interested in the set of all motifs
		PssmMotif::set_t all_motifs_set;
		std::copy(pssm_motifs.begin(), pssm_motifs.end(), std::inserter(all_motifs_set, all_motifs_set.begin()));
		motif_sets.insert(all_motifs_set);

		//are we interested in every subset of all the motifs with size N-1?
		if (leave_one_out)
		{
			for (PssmMotif::vec_t::const_iterator m = pssm_motifs.begin();
				pssm_motifs.end() != m;
				++m)
			{
				PssmMotif::set_t one_left_out = all_motifs_set;
				one_left_out.erase(*m);
				motif_sets.insert(one_left_out);
			}
		}
	}



	void add_score_to_ranks(const PssmMotif::score_map_t & score_map, bool is_sequence_score, const std::string & name)
	{
		//for each set of motifs we're interested in
		for (motif_set_set_t::const_iterator s = motif_sets.begin();
			motif_sets.end() != s;
			++s)
		{
			//build the score map for this set of motifs
			PssmMotif::score_map_t scores;
			for (PssmMotif::score_map_t::const_iterator i = score_map.begin();
				score_map.end() != i;
				++i)
			{
				if (s->find(i->first) != s->end())
				{
					scores.insert(*i);
				}
			}

			//what is the score?
			const double score = PssmMotif::get_score(scores);

			//which rank map do we want to put the score in?
			rank_map_t & r = is_sequence_score ? rank_map[*s].get<0>() : rank_map[*s].get<1>();

			//insert the score
			r.insert(
				rank_map_t::value_type(
					BIO_NS::float_t( score ),
					name));
		}
	}

	void display_best_ranked() const
	{
		//for each set of motifs we're interested in
		for (motif_set_set_t::const_iterator s = motif_sets.begin();
			motif_sets.end() != s;
			++s)
		{
			std::cout << "Displaying results for motif set:\n";
			for (PssmMotif::set_t::const_iterator i = s->begin();
				s->end() != i;
				++i)
			{
				std::cout << motif_descriptions_map.find(*i)->second << "\n";
			}

			//the ranking for this motif set
			const motif_ranking_t & motif_ranking = rank_map.find(*s)->second;

			//display the best sequences
			unsigned num_displayed = 0;
			cout << "Best sequence groups:\n";
			for (rank_map_t::const_reverse_iterator r = motif_ranking.get<0>().rbegin();
				motif_ranking.get<0>().rend() != r && num_displayed != num_to_display;
				++r)
			{
				cout << r->first << '\t' << r->second << '\n';
				++num_displayed;
			}
			cout << "\n";

			//display the best remos
			cout << "Best ranges:\n";
			num_displayed = 0;
			for (rank_map_t::const_reverse_iterator r = motif_ranking.get<1>().rbegin();
				motif_ranking.get<1>().rend() != r && num_displayed != num_to_display;
				++r)
			{
				cout << r->first << '\t' << r->second << '\n';
				++num_displayed;
			}
			cout << "\n";
		}
	}

	int task()
	{
		if (motif_descriptions.empty())
		{
			throw std::logic_error( "No motifs specified" );
		}

		parse_motif_descriptions();

		build_interesting_motif_sets();

		deserialise_analysis();

		cout << "\nLooking for the following motifs:\n";
		copy(motif_descriptions.begin(), motif_descriptions.end(), ostream_iterator< std::string >(cout, "\n"));

		//load biobase
		BiobaseDb::singleton();

		visit_remo_analysis();

		display_best_ranked();

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return PssmMotifRankerApp().main(argc, argv);
}

