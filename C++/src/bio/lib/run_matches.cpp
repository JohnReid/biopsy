/* Copyright John Reid 2007-2010
*/

#include "bio-pch.h"


#include "bio/defs.h"



#include "bio/open_file.h"
#include "bio/sequence.h"
#include "bio/biobase_filter.h"
#include "bio/svg_match.h"
#include "bio/environment.h"
#include "bio/pssm_match.h"
#include "bio/run_match.h"
#include "bio/biobase_data_traits.h"
#include "bio/svg_match.h"
#include "bio/bifa_analysis.h"

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <xercesc/util/PlatformUtils.hpp>
XERCES_CPP_NAMESPACE_USE

#include <string>
#include <set>
using namespace std;

BIO_NS_START


void
find_factors_for_matches(
	match_result_vec_t & hits,
	factor_scores_map_t & factor_scores)
{
	for (match_result_vec_t::iterator h = hits.begin();
		hits.end() != h;
		++h)
	{
		const BiobaseTablePssmEntry * entry = BiobaseDb::singleton().get_pssm_entry(h->link);
		const FactorLinkList & factors = entry->get_factors();

		//we give each factor hit by this matrix a fraction of the hits total value (1.0)
		const float_t score_per_factor =
			1.0 == h->result.score
				? std::numeric_limits<float_t>::max()
				: float_t(1.0) / float_t(factors.size()) / (float_t(1.0) - float_t(h->result.score));

		for (FactorLinkList::const_iterator f = factors.begin();
			factors.end() != f;
			++f)
		{
			FactorInfo & factor_info = factor_scores[(*f)->link];
			factor_info.score += score_per_factor;
			factor_info.hits.insert(h->number);
		}
	}
}

BuildSvgArgs::BuildSvgArgs(
	const std::string & file,
	const std::string & title,
	float_t max_threshold,
	float_t min_threshold,
	unsigned max_num_factors,
	bool show_labels,
	bool open_file,
	const std::string & notes )
	: file( file )
	, title( title )
	, max_threshold( max_threshold )
	, min_threshold( min_threshold )
	, max_num_factors( max_num_factors )
	, show_labels( show_labels )
	, open_file( open_file )
	, notes( notes )
{ }

void
BuildSvgArgs::add_options(boost::program_options::options_description & options)
{
	options.add_options()
		("min_display", po::value(&min_threshold)->default_value(0.0f), "min threshold for hits displayed")
		("max_display", po::value(&max_threshold)->default_value(1.0f), "max threshold for hits displayed")
		("svg_file", po::value(&file)->default_value("bifa.svg"), "SVG file to write to")
		("svg_title", po::value(&title)->default_value("BiFa\nAnalysis"), "title for SVG")
		("max_num_factors", po::value(&max_num_factors)->default_value(16), "maximum # factors to show")
		("show_svg_labels", po::value(&show_labels)->default_value(false), "show SVG labels")
		("open_svg", po::bool_switch(&open_file), "open SVG file")
		("notes", po::value(&notes)->default_value(""), "notes to add to SVG")
		;
}


/** Run the algorithm with the arguments. */
void
BuildSvgArgs::run_build_svg(
	const seq_t & seq,
	match_result_vec_t & results)
{
	build_svg(
		boost::filesystem::path(file),
		title,
		seq,
		min_threshold,
		results,
		max_num_factors,
		show_labels,
		open_file,
		0,
		notes );
}

std::ostream &
operator<<(std::ostream & os, const BuildSvgArgs & args)
{
	os
		<< "Min threshold to display: " << args.min_threshold << "\n"
		<< "Max threshold to display: " << args.max_threshold << "\n"
		<< "SVG file: " << args.file << "\n"
		<< "SVG title: " << args.title << "\n"
		<< "Max # factors: " << args.max_num_factors << "\n"
		<< "Show SVG labels: " << (args.show_labels ? "yes" : "no") << "\n"
		<< "Open SVG file: " << (args.open_file ? "yes" : "no") << "\n"
		<< "Notes: " << args.notes << "\n"
		;

	return os;
}

namespace detail {
struct hit_is_same
	: std::binary_function< MatchResults, MatchResults, bool >
{
	bool operator()( const MatchResults & lhs, const MatchResults & rhs ) const
	{
		return
			lhs.link == rhs.link
			&&
			lhs.result.position == rhs.result.position
			&&
			lhs.result.complement == rhs.result.complement;
	}
};
}

void
build_svg(
	const fs::path & file,
	const string & title,
	const seq_t & seq,
	float_t min_threshold,
	match_result_vec_t & results,
	size_t max_num_factors,
	bool show_labels,
	bool open,
	match_result_vec_t * max_chain,
	const std::string & notes,
	float max_threshold )
{
	XMLPlatformUtils::Initialize();

	//find the maximum score
	for (match_result_vec_t::const_iterator i = results.begin();
		results.end() != i;
		++i)
	{
		max_threshold = std::max( max_threshold, i->result.score );
	}

	SvgDomDocument doc(seq.size(), min_threshold, max_threshold, title, seq, show_labels);

	//for each hit result
	size_t idx = 0;
	unsigned num_added = 0;
	for (match_result_vec_t::iterator i = results.begin();
		results.end() != i;
		++i, ++idx)
	{
		i->number = idx;

		BiobaseTablePssmEntry * pssm_entry = BiobaseDb::singleton().get_pssm_entry(i->link);

		assert((size_t) i->result.position < seq.size());

		if (min_threshold <= i->result.score)
		{
			const bool is_in_max_chain =
				0 != max_chain
				&&
				max_chain->end() != std::find_if( max_chain->begin(), max_chain->end(), std::bind1st( detail::hit_is_same(), *i ) );

			doc.add_result(
				i->result,
				*pssm_entry,
				idx,
				is_in_max_chain );
			++num_added;
		}
	}
	if (0 == num_added)
	{
		throw std::logic_error( "No hits above threshold" );
	}

	factor_scores_map_t factors;
	find_factors_for_matches(results, factors);
	//cout << factors.size() << " factors associated with these matches" << std::endl;

	//for each factor - in score order
	while (! factors.empty() && 0 != max_num_factors--)
	{
		//find the factor with the highest score
		factor_scores_map_t::iterator best = factors.end();
		for (factor_scores_map_t::iterator f = factors.begin();
			factors.end() != f;
			++f)
		{
			if (factors.end() == best || f->second.score > best->second.score)
			{
				best = f;
			}
		}
		assert(best != factors.end());

		Factor * factor = BiobaseDb::singleton().get_entry<FACTOR_DATA>(best->first);
		if (0 != factor)
		{
			doc.add_factor(factor, best->second.hits);
		}

		factors.erase(best);
	}

	if( notes.size() )
	{
		doc.add_notes( notes );
	}

	add_tooltip_support(doc.doc, doc.doc_root);

	dom_print(doc.doc, file._BOOST_FS_NATIVE().c_str());

	//copy the script to the same directory
	try
	{
		fs::path script(BioEnvironment::singleton().get_svg_script_file());
		fs::path script_dest(file.branch_path());
		script_dest /= script.leaf();

		//only copy if destination older
		if (! fs::exists(script_dest) || fs::last_write_time(script) > fs::last_write_time(script_dest))
		{
			//remove first if exists
			if (fs::exists(script_dest))
			{
				fs::remove(script_dest);
			}

			fs::copy_file(script, script_dest);
		}
	}
	catch (const std::exception & ex)
	{
		cerr << ex.what() << endl;
	}

	XMLPlatformUtils::Terminate();

	if (open)
	{
		open_file(file);
	}

}


void
pssm_match(
	const seq_t & match_seq,
	float_t threshold,
	ScoreAlgorithm algorithm,
	const fs::path & file,
	const std::string & title,
	bool show_labels)
{
#ifdef TEST
	set<string> species_acronyms;
	for (Site::map_t::const_iterator i = BiobaseDb::singleton().get_sites().begin();
		BiobaseDb::singleton().get_sites().end() != i;
		++i)
	{
		species_acronyms.insert(i->second->id.species_group);
	}
	copy(species_acronyms.begin(), species_acronyms.end(), ostream_iterator<string>(cout, ", "));
	cout << endl;
#endif

	BiobasePssmFilter filter; //all consensus sequences & vertebrates

	//score the sites and matrices
	match_result_vec_t results;
	pssm_match(
		match_seq.begin(),
		match_seq.end(),
		threshold,
		algorithm,
		false,
		filter,
		filter,
		inserter(results, results.begin()));

#ifdef VERBOSE_CHECKING
	cout << results.size() << " matches over threshold" << endl;
#endif

	cout << "Building SVG" << endl;
	build_svg(
		file,
		title,
		match_seq,
		threshold,
		results,
		16,
		show_labels,
		true );
}



BIO_NS_END
