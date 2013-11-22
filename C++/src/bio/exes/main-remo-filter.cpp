/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo_analysis.h"
#include "bio/tss_data.h"
#include "bio/environment.h"
USING_BIO_NS

#include <boost/progress.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;





struct FilterReMosApp : Application
{
	std::string remo_extraction_filename;
	std::string output_filename;
	bool filter_on_clones;
	bool allow_tigr_clones;
	bool remove_old_analyses;
	bool only_upstream;
	bool remove_exons;
	double allowed_exon_overlap;
	unsigned min_num_species;
	ReMoExtraction::ptr_t remo;
	std::string sequence;

	FilterReMosApp()
	{
		get_options().add_options()
			("remo,r", po::value(&remo_extraction_filename)->default_value("remo_extraction.bin"), "remo extraction file")
			("output,o", po::value(&output_filename)->default_value("remo_space.bin"), "output file")
			("filter_on_clones,c", po::value(&filter_on_clones)->default_value(false), "filter remos without RIKEN clones")
			("allow_tigr_clones,t", po::value(&allow_tigr_clones)->default_value(false), "allow TIGR clones")
			("only_upstream,u", po::value(&only_upstream)->default_value(false), "filter upstream")
			("filter_exons,e", po::value(&remove_exons)->default_value(false), "filter exons")
			("allowed_exon_overlap,a", po::value(&allowed_exon_overlap)->default_value(0.1), "allowed exon overlap")
			("filter_old_analyses,d", po::value(&remove_old_analyses)->default_value(true), "filter old dbs")
			("min_num_species,m", po::value(&min_num_species)->default_value(0), "filter remos in fewer species")
			("sequence,s", po::value(&sequence)->default_value(""), "regex to match sequence name")
			;
	}

	int task()
	{
		cout << "Reading remos from " << remo_extraction_filename << "\n";
		cout << (filter_on_clones ? "Filtering" : "Not filtering") << " remos without TSS clone data\n";
		cout << (allow_tigr_clones ? "Using" : "Ignoring") << " TIGR clone data\n";
		cout << (only_upstream ? "Filtering" : "Not filtering") << " non-upstream remos\n";
		if (remove_exons)
		{
			cout << "Filtering remos that overlap exons more than " << int(allowed_exon_overlap * 100) << "%\n";
		}
		else
		{
			cout << "Not filtering remos that overlap exons\n";
		}
		cout << (remove_old_analyses ? "Filtering" : "Not filtering") << " remos using older db versions\n";
		if (min_num_species > 0)
		{
			cout << "Filtering remos in fewer than " << min_num_species << " species\n";
		}
		else
		{
			cout << "Not filtering remos based on number of species\n";
		}
		if ("" != sequence)
		{
			cout << "Filtering remos that don't match regex: " << sequence << "\n";
		}
		else
		{
			cout << "Not filtering remos based on name\n";
		}

		//deserialise the binary remo archive
		{
			fs::path
				remo_extraction_archive(
					remo_extraction_filename
				);
			cout << "Deserialising remo extraction from \"" << remo_extraction_archive._BOOST_FS_NATIVE() << "\"\n";
			boost::progress_timer timer;
			remo = ReMoExtraction::deserialise(remo_extraction_archive);
		}

		//remove initial empty sequence groups
		cout << "*********** Removing empty sequence groups ***********\n";
		AlwaysFalsePred always_false_pred;
		erase_if(always_false_pred);

		if ("" != sequence)
		{
			cout << "*********** Removing remos that don't match regex ***********\n";
			regex re(sequence);
			SequenceRegexPred pred(re);
			erase_if(pred);
		}

		//do we want to only use those from the most recent dbs?
		if (remove_old_analyses)
		{
			cout << "*********** Removing remos from comparisons using out of date genomes **********\n";
			OldDbPred pred(remo);
			erase_if(pred);
		}

		//do we want to remove remos that intersect exons?
		if (remove_exons)
		{
			cout << "*********** Removing remos that intersect exons **********\n";
			ExonPred pred(allowed_exon_overlap);
			erase_if(pred);
		}

		//do we want to only use upstream remos?
		if (only_upstream)
		{
			cout << "*********** Removing non-upstream remos **********\n";
			NotUpstreamPred pred;
			erase_if(pred);
		}

		//do we want to only use those with clones?
		if (filter_on_clones)
		{
			cout << "*********** Removing remos without clone TSS data ***********\n";
			cout
				<< (allow_tigr_clones
					? "Looking for RIKEN & TIGR clones\n"
					: "Only looking for RIKEN clones\n");

			ClonePred pred(allow_tigr_clones);
			erase_if(pred);
			cout
				<< "Could not find TSS data for " << pred.num_remos_without_tss_data << " remos in "
				<< pred.genes_without_tss_data.size() << " genes\n"
				<< "Found " << pred.num_with_riken_clone << " with RIKEN clones\n"
				<< "Found " << pred.num_with_tigr_clone << " with TIGR clones\n\n";
		}

		//do we want to only use those remos with a minimum # of sequences?
		if (min_num_species > 1)
		{
			cout << "*********** Removing remos in fewer than " << min_num_species << " species ***********\n";
			MinNumSpeciesPred pred(min_num_species);
			erase_if(pred);
		}

		unsigned total_num_groups = 0;
		unsigned total_num_remos = 0;
		//count the number of groups and remos
		for (ReMoSequenceGroup::list_t::iterator sg = remo->sequence_groups.begin();
			remo->sequence_groups.end() != sg;
			++sg, ++total_num_groups)
		{
			for (ReMoBundle::map_t::iterator rb = sg->get()->remo_bundles.begin();
				sg->get()->remo_bundles.end() != rb;
				++rb, ++total_num_remos)
			{
			}
		}
		cout << "\nLeft with " << total_num_groups << " groups containing " << total_num_remos << " remos\n\n";

		//serialise the reduced remo map
		{
			cout << "Serialising remo extraction part to \"" << output_filename << "\"\n";
			std::ofstream stream(output_filename.c_str(), std::ios::binary);
			boost::archive::binary_oarchive(stream) << const_cast<const ReMoExtraction &>(*remo);
		}

		return 0;
	}

	/** Erase all those remos for which the predicate returns true. Also erase empty sequence groups. */
	template <typename Pred>
	void erase_if(Pred & pred)
	{
		unsigned num_remos_erased = 0;
		unsigned total_num_remos = 0;
		unsigned num_groups_erased = 0;
		unsigned total_num_groups = 0;
		unsigned num_centre_sequences_erased = 0;
		std::set< std::string > genes_before;
		std::set< std::string > genes_after;

		for (ReMoSequenceGroup::list_t::iterator sg = remo->sequence_groups.begin();
			remo->sequence_groups.end() != sg;
			++total_num_groups)
		{
			std::string gene_id = get_gene_id(sg->get()->get_centre_sequence().id);
			genes_before.insert(gene_id);

			for (ReMoBundle::map_t::iterator rb = sg->get()->remo_bundles.begin();
				sg->get()->remo_bundles.end() != rb;
				++total_num_remos)
			{
				//if the predicate is true
				if (pred(*sg, *rb))
				{
					sg->get()->remo_bundles.erase(rb++);
					++num_remos_erased;
				}
				else
				{
					//check to see if there is a centre sequence
					if (rb->second->remos.end() == rb->second->remos.find(rb->second->centre_sequence))
					{
						//there isn't
						sg->get()->remo_bundles.erase(rb++);
						++num_centre_sequences_erased;
						++num_remos_erased;
					}
					else
					{
						++rb;
					}
				}
			}

			if (sg->get()->remo_bundles.empty())
			{
				remo->sequence_groups.erase(sg++);
				++num_groups_erased;
			}
			else
			{
				genes_after.insert(gene_id);
				++sg;
			}
		}
		cout << "Missing " << num_centre_sequences_erased << " centre sequences\n";
		cout << "Started with " << genes_before.size() << " genes, left with " << genes_after.size() << "\n";
		cout << "Discarded " << num_remos_erased << " remos from total of " << total_num_remos << "\n";
		cout << "Discarded " << num_groups_erased << " empty sequence groups from total of " << total_num_groups << "\n\n\n";
	}

	/** Used to remove initial empty sequence groups. */
	struct AlwaysFalsePred
	{
		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			return false;
		}
	};

	/** True if remo is not upstream. */
	struct NotUpstreamPred
	{
		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			//cout << group->sequences.begin()->get()->location << "\n";
			return group->sequences.begin()->get()->location != REMO_LOC_UPSTREAM;
		}
	};

	/** True if regex doesn't match sequence name. */
	struct SequenceRegexPred
	{
		regex re;

		SequenceRegexPred(regex re) : re(re) { }

		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			smatch what;
			bool matched = false;
			//check each of the sequences in the group to see if the id matches
			for (ReMoSequence::list_t::const_iterator s = group->sequences.begin();
				group->sequences.end() != s && ! matched;
				++s)
			{
				matched = regex_search(s->get()->id, what, re);
			}

			return ! matched;
		}
	};

	/** True if remo is intersects an exon. */
	struct MinNumSpeciesPred
	{
		unsigned min_num_species;

		MinNumSpeciesPred(unsigned min_num_species)
			: min_num_species(min_num_species)
		{
		}

		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			return remo.second->remos.size() < min_num_species;
		}
	};

	/** True if remo is intersects an exon. */
	struct ExonPred
	{
		double allowed_overlap;

		ExonPred(double allowed_overlap)
			: allowed_overlap(allowed_overlap)
		{
		}

		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			//for each species comprising this bundle
			for (ReMo::map_t::iterator r = remo.second->remos.begin();
				remo.second->remos.end() != r;
				)
			{
				//get the sequence for it
				ReMoSequence::ptr_t sequence = group->get_sequence_for(r->first);

				//for each part of the remo in this species
				for (ReMo::list_t::iterator p = r->second.begin();
					r->second.end() != p;
					)
				{
					//check the exons to see if they overlap it
					ReMoExon::list_t::iterator e = sequence->exons.begin();
					for ( ;
						sequence->exons.end() != e;
						++e)
					{
						if (p->get()->range.overlap(e->range) > allowed_overlap)
						{
							//they do overlap - so remove this part of this remo in this species
							r->second.erase(p++);
							break; //out of exon iteration
						}
					}

					//did we delete the part?
					if (sequence->exons.end() == e)
					{
						//no
						++p;
					}
				}

				//did we remove all of the parts in this species?
				if (r->second.empty())
				{
					//yes - so remove the species entry in the map
					remo.second->remos.erase(r++);
				}
				else
				{
					//no - carry on
					++r;
				}
			}

			//return true (i.e. to erase) if there is only one species or less left
			return remo.second->remos.size() < 2;
		}
	};

	struct OldDbPred
	{
		typedef map< string, string > gene_versions_map_t;
		gene_versions_map_t gene_versions_map;

		OldDbPred(ReMoExtraction::ptr_t remo)
		{
			//iterate through the extraction
			{
				cout << "Parsing remo ids\n";
				boost::progress_timer timer;

				typedef map< string, set< string > > species_versions_map_t;
				species_versions_map_t species_versions_map;

				typedef map< string, unsigned > species_versions_counts_t;
				map<string, unsigned> species_versions_counts;

				unsigned unparsed_count = 0;

				for (ReMoSequenceGroup::list_t::const_iterator sg = remo->sequence_groups.begin();
					remo->sequence_groups.end() != sg;
					++sg)
				{
					for (ReMoBundle::map_t::const_iterator rb = sg->get()->remo_bundles.begin();
						sg->get()->remo_bundles.end() != rb;
						++rb)
					{
						ReMoSequenceId id;
						if (parse_sequence_id(rb->second->centre_sequence, id))
						{
							map<string, unsigned>::iterator c = species_versions_counts.find(id.species);
							if (species_versions_counts.end() == c)
							{
								c = species_versions_counts.insert(map<string, unsigned>::value_type(id.species, 0)).first;
							}
							c->second++;
							
							typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
							boost::char_separator<char> sep("_");
							tokenizer tokens(id.species, sep);
							vector<string> words;
							copy(tokens.begin(), tokens.end(), inserter(words, words.begin()));
							if (5 == words.size() && "core" == words[2])
							{
								const string species = words[0] + " " + words[1];
								const string version = words[3] + "_" + words[4];
								species_versions_map[species].insert(version);
							}
		
							const string gene = id.gene_id + " " + id.transcript_id;
							if ("" == gene_versions_map[gene] || gene_versions_map[gene] < id.species)
							{
								gene_versions_map[gene] = id.species;
							}
						}
						else
						{
							unparsed_count++;
						}
					}
				}

				cout << "Failed to parse " << unparsed_count << " species ids\n";

				cout << "Found these species and versions\n";
				for (species_versions_map_t::const_iterator s = species_versions_map.begin();
					species_versions_map.end() != s;
					++s)
				{
					cout << s->first << ": ";
					for (set<string>::const_iterator v = s->second.begin();
						s->second.end() != v;
						++v)
					{
						cout << *v << ", ";
					}
					cout << "\n";
				}
		
				cout << "this number of times\n";
				for (species_versions_counts_t::const_iterator c = species_versions_counts.begin();
					species_versions_counts.end() != c;
					++c)
				{
					cout << c->first << ": " << c->second << "\n";
				}
			}
		}

		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			bool to_remove = true;

			ReMoSequenceId id;
			if (parse_sequence_id(remo.second->centre_sequence, id))
			{
                const string gene = id.gene_id + " " + id.transcript_id;
				if (gene_versions_map[gene] == id.species)
				{
					//don't remove it if it is the most recent database version
					to_remove = false;
				}
			}

			return to_remove;
		}
	};

	struct ClonePred
	{
		TSS::map_t tss_map;
		unsigned num_remos_without_tss_data;
		unsigned num_with_riken_clone;
		unsigned num_with_tigr_clone;
		std::set< std::string > genes_without_tss_data;
		bool allow_tigr_clones;

		ClonePred(bool allow_tigr_clones = false)
			: num_remos_without_tss_data(0)
			, num_with_riken_clone(0)
			, num_with_tigr_clone(0)
			, allow_tigr_clones(allow_tigr_clones)
		{
			//load the TSS data
			{
				cout
					<< "Loading TSS data from \""
					<< BioEnvironment::singleton().get_tss_file()
					<< "\" and \""
					<< BioEnvironment::singleton().get_tss_clones_file() << "\"\n";
				//boost::progress_timer timer;

				fs::path tss_file(BioEnvironment::singleton().get_tss_file());
				fs::path clones_file(BioEnvironment::singleton().get_tss_clones_file());
				TSS::parse_files(tss_file, clones_file, tss_map);
				cout << "Loaded TSS data for " << tss_map.size() << " genes\n";
			}
		}

		bool operator()(ReMoSequenceGroup::list_t::value_type & group, ReMoBundle::map_t::value_type & remo)
		{
			ReMoSequenceId id;
			if (parse_sequence_id(remo.second->centre_sequence, id))
			{
				//find the TSS data for this gene
				TSS::map_t::const_iterator t = tss_map.find(id.gene_id);
				// cout << id.gene_id << " - " << tss_map.begin()->first << "\n";
				if (tss_map.end() == t)
				{
					++num_remos_without_tss_data;
					genes_without_tss_data.insert(id.gene_id);
				}
				else
				{
					const TSS & tss = t->second;

					//for each clone we have, check if it is riken clone of right transcript
					bool found_clone = false;
					for (Clone::vec_t::const_iterator c = tss.clones.begin();
						tss.clones.end() != c;
						++c)
					{
						// cout << c->transcript << " - " << id.transcript_id << "\n";
						if (c->transcript == id.transcript_id)
						{
							//check for RIKEN clone first
							if (RIKEN_CLONE == c->type)
							{
								found_clone = true;
								++num_with_riken_clone;
								break;
							}
							//now check for TIGR clones
							else if (allow_tigr_clones && TIGR_CLONE == c->type)
							{
								//if the TIGR clone is
								found_clone = true;
								++num_with_tigr_clone;
								break;
							}
						}
					}
					if (found_clone)
					{
						return false;
					}
				}
			}

			return true;
		}
	};

};

int
main(int argc, char * argv[])
{
	return FilterReMosApp().main(argc, argv);
}

