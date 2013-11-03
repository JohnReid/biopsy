/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/remo_analysis.h"
#include "bio/binding_hit.h"
#include "bio/binding_model.h"

#include <boost/regex.hpp>
#include <boost/progress.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/version.hpp>
namespace po = boost::program_options;

#include <fstream>
#include <iostream>
using namespace std;

BIO_NS_START

void
ReMoAnalysis::convert(map_t & remo_map, BiFaAnalysis::map_t & bifa_map)
{
	for (map_t::const_iterator s = remo_map.begin();
		remo_map.end() != s;
		++s)
	{
		for (location_map_t::const_iterator l = s->second.begin();
			s->second.end() != l;
			++l)
		{
			for (sequence_map_t::const_iterator r = l->second.begin();
				l->second.end() != r;
				++r)
			{
				std::stringstream key;
				key << s->first << " " << l->first << " " << r->first;

				bifa_map[key.str()].reset(new BiFaAnalysis);
				bifa_map[key.str()]->sequence = r->second->sequence;
				std::copy(
					r->second->results.begin(),
					r->second->results.end(),
					std::inserter( bifa_map[key.str()]->results, bifa_map[key.str()]->results.begin() ) );
			}
		}
	}
}

void
AnalysisVisitor::add_analysis_options(boost::program_options::options_description & options)
{
	options.add_options()
		("analysis,a", po::value(&analysis_filename), "input analysis file")
		("bifa,b", po::value(&is_bifa_analysis)->default_value(false), "is bifa analysis?")
		("sequence_name_regex,s", po::value(&sequence_name_regex), "regex for sequence names")
		
		;
}

bool
AnalysisVisitor::visit_sequence_group(const std::string & seq_group_name)
{
	return true;
}

void
AnalysisVisitor::leave_sequence_group(const std::string & seq_group_name)
{
}

void
AnalysisVisitor::visit_remo(
	const std::string & seq_group_name,
	ReMoLocation location,
	const ReMoRange & range,
	const std::string & remo_name,
	bifa_hits_t & results,
	const seq_t & sequence)
{
}

void
AnalysisVisitor::deserialise_analysis()
{
	if ("" == analysis_filename)
	{
		throw std::logic_error( "No analysis file specified" );
	}

	cout << "Deserialising " << (is_bifa_analysis ? "bifa" : "extraction") << " analysis from \"" << analysis_filename << "\"\n";
	std::ifstream stream(analysis_filename.c_str(), std::ios::binary);

	if (is_bifa_analysis)
	{
		boost::archive::binary_iarchive(stream) >> bifa_analysis;
	}
	else
	{
		boost::archive::binary_iarchive(stream) >> remo_analysis;
	}
}

void
AnalysisVisitor::serialise_analysis(const std::string & filename) const
{
	if ("" == filename)
	{
		throw std::logic_error( "No file specified" );
	}

	cout << "Serialising analysis to \"" << filename << "\"\n";
	std::ofstream stream(filename.c_str(), std::ios::binary);

	if (is_bifa_analysis)
	{
		boost::archive::binary_oarchive(stream) << bifa_analysis;
	}
	else
	{
		boost::archive::binary_oarchive(stream) << remo_analysis;
	}
}


void
AnalysisVisitor::visit_remo_analysis(bool show_progress)
{
	num_sequences = 0;
	num_remos = 0;
	num_bases = 0;
	num_hits = 0;
	expected_num_hits = 0.0;

	boost::smatch what;
	boost::regex seq_re;
	if ("" == sequence_name_regex)
	{
		std::cout << "Visiting remos\n";
		seq_re = boost::regex(".");
	}
	else
	{
		std::cout << "Visiting remos that match regex: \"" << sequence_name_regex << "\"\n";
		seq_re = boost::regex(sequence_name_regex);
	}

	if (is_bifa_analysis)
	{
		boost::scoped_ptr< boost::progress_display > sp(show_progress ? new boost::progress_display( bifa_analysis.size()) : 0);

		for (BiFaAnalysis::map_t::iterator s = bifa_analysis.begin();
			bifa_analysis.end() != s;
			)
		{
			0 == sp || ++(*sp);

			//does it match the regex and do we want to visit this sequence group?
			if ((! boost::regex_search(s->first, what, seq_re)) || (! visit_sequence_group(s->first)))
			{
				//no
				bifa_analysis.erase(s++);
				continue;
			}

			++num_sequences;
			++num_remos;
			num_bases += s->second->sequence.size();
			num_hits += s->second->results.size();
			expected_num_hits += get_expected_num_hits( s->second->results.begin(), s->second->results.end() );

			visit_remo(
				s->first,
				REMO_LOC_UNDEFINED,
				ReMoRange(0, s->second->sequence.length()),
				s->first,
				s->second->results,
				s->second->sequence);

			leave_sequence_group(s->first);

			++s;
		}
	}
	else
	{
		boost::scoped_ptr< boost::progress_display > sp(show_progress ? new boost::progress_display( remo_analysis.size()) : 0);
		for (ReMoAnalysis::map_t::iterator s = remo_analysis.begin();
			remo_analysis.end() != s;
			)
		{
			0 == sp || ++(*sp);

			//does it match the regex and do we want to visit this sequence group?
			if ((! boost::regex_search(s->first, what, seq_re)) || (! visit_sequence_group(s->first)))
			{
				//no
				remo_analysis.erase(s++);
				continue;
			}

			++num_sequences;

			//for each location
			for (ReMoAnalysis::location_map_t::const_iterator l = s->second.begin();
				s->second.end() != l;
				++l)
			{
				//for each remo
				for (ReMoAnalysis::sequence_map_t::const_iterator r = l->second.begin();
					l->second.end() != r;
					++r)
				{
					visit_remo(
						s->first,
						l->first,
						r->first,
						BIO_MAKE_STRING(s->first << " " << l->first << " " << r->first),
						r->second->results,
						r->second->sequence);

					++num_remos;
					num_hits += r->second->results.size();
					num_bases += r->second->sequence.size();
					expected_num_hits += get_expected_num_hits( r->second->results.begin(), r->second->results.end() );

				}
			}

			leave_sequence_group(s->first);

			++s;
		}
	}

	std::cout
		<< "Visited " << num_sequences
		<< " sequence groups with " << num_remos
		<< " remos, " << num_bases
		<< " bases, " << num_hits
		<< " putative hits and " << expected_num_hits
		<< " expected hits\n";
}



BIO_NS_END

#ifndef NDEBUG
namespace boost { namespace serialization {
/// work-around for undefined symbol in debug mode using boost SVN revision 64053.
//template <>
//const unsigned int
//version< boost::multi_index::detail::serialization_version< BIO_NS::BindingHit< BIO_NS::BindingModel > > >::value
//	= version< boost::multi_index::detail::serialization_version< BIO_NS::BindingHit< BIO_NS::BindingModel > > >::type::value
//	;
/// Explicit instantiation
template
struct version< boost::multi_index::detail::serialization_version< BIO_NS::BindingHit< BIO_NS::BindingModel > > >;

/// Explicit instantiation
//template
//const unsigned int version< boost::multi_index::detail::serialization_version< BIO_NS::BindingHit< BIO_NS::BindingModel > > >::value;
} }
#endif //NDEBUG

//boost::serialization::version<boost::multi_index::detail::serialization_version<bio::BindingHit<bio::BindingModel> > >::value


