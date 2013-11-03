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



struct PssmMotifFinderApp : Application, AnalysisVisitor
{
	typedef std::vector< std::string > string_vec_t;

	std::string output_analysis_filename;
	string_vec_t motif_descriptions;
	bool show_matches;
	PssmMotif::vec_t pssm_motifs;

	PssmMotifFinderApp()
	{
		get_options().add_options()
			("output,o", po::value(&output_analysis_filename), "output analysis file")
			("show_matches,m", po::value(&show_matches)->default_value(false), "show motif matches")
			("motif", po::value(&motif_descriptions), "motif description")
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
		}
	}

	void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		bifa_hits_t & hits,
		const seq_t & sequence)
	{
		match_result_vec_t results;
		bifa_hits_2_match_results( hits, results );
		if( results.empty() )
		{
			return;
		}

		//cout << "Searching analysis of " << remo_name << "\n";

		typedef std::set<MatchResults> hit_set_t;
		hit_set_t hits_in_motifs_set;

		//for each motif
		BIO_NS::float_t min_score = 1.0f;
		for (PssmMotif::vec_t::const_iterator m = pssm_motifs.begin();
			pssm_motifs.end() != m;
			++m)
		{
			//look for the motifs
			PssmMotif::HitVec hits;
			m->get()->find_in( results, hits );

			//for each hit
			for (PssmMotif::HitVec::const_iterator hit = hits.begin();
				hits.end() != hit;
				++hit)
			{
				for (PssmMotif::HitElement::vec_t::const_iterator h = hit->begin();
					hit->end() != h;
					++h)
				{
					hits_in_motifs_set.insert( *( h->match_result ) );
					min_score = std::min( h->match_result->result.score, min_score);
				}
				if (show_matches)
				{
					cout << *hit << "\n";
				}
			}
		}

		if( ! hits_in_motifs_set.empty() )
		{
			std::cout << remo_name << std::endl;
		}
	}




	int task()
	{
		if (motif_descriptions.empty())
		{
			throw std::logic_error( "No motifs specified" );
		}

		parse_motif_descriptions();

		deserialise_analysis();

		cout << "\nLooking for the following motifs:\n";
		copy(motif_descriptions.begin(), motif_descriptions.end(), ostream_iterator< std::string >(cout, "\n"));

		//load biobase
		BiobaseDb::singleton();

		visit_remo_analysis( false );

		//do we want to write the output?
		if ("" != output_analysis_filename)
		{
			//serialise
			serialise_analysis(output_analysis_filename);
		}

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return PssmMotifFinderApp().main(argc, argv);
}

