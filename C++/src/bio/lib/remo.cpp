/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/remo.h"
#include "bio/lexer.h"

#include "ReMoParser.hpp"
#include "ReMoLexer.hpp"
#include "ReMoNumberLexer.hpp"
#include "ReMoStringLexer.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
namespace fs = boost::filesystem;

#include <algorithm>



BIO_NS_START

class ReMoSharedLexer
{
public:
	ReMoSharedLexer(std::istream & is)
		: remoLexer(is)
		, numberLexer(remoLexer.getInputState())
		, stringLexer(remoLexer.getInputState())
	{
		selector.addInputStream(&remoLexer, "remo");
		selector.select("remo"); //start state

		selector.addInputStream(&numberLexer, "number");
		selector.addInputStream(&stringLexer, "string");
	}

	antlr::TokenStreamSelector & get_selector() { return selector; }

protected:
	ReMoLexer						remoLexer;
	ReMoNumberLexer					numberLexer;
	ReMoStringLexer					stringLexer;
	antlr::TokenStreamSelector		selector;
};


ReMoExtraction::ptr_t
parse_remo_extraction(boost::filesystem::path file)
{
	ReMoExtraction::ptr_t result;

	fs::ifstream is(file);
	if (! is)
	{

		throw BIO_MAKE_STRING("Could not open file: " << file._BOOST_FS_NATIVE());
	
	}
	else
	{
	
		ReMoSharedLexer lexer(is);

		ReMoParser parser(lexer.get_selector());
		parser.sps.selector = &lexer.get_selector();

		result = parser.remo_extraction();
	}

	return result;
}

ReMoSpecies parse_remo_species(const std::string & species_string)
{
	return REMO_SPECIES_CHICK;
}

ReMoSubSequence
ReMoBundle::get_centre_sequence_id() const
{
	if (remos.empty())
	{
		throw std::logic_error( "No sequences in bundle" );
	}

	ReMo::map_t::const_iterator i = remos.find(centre_sequence);
	BOOST_ASSERT(remos.end() != i);
	BOOST_ASSERT(! i->second.empty());
	return ReMoSubSequence(centre_sequence, i->second.begin()->get()->range);
}

ReMoBundle::id_set_t
ReMoBundle::get_sequence_ids() const
{
	id_set_t result;

	for (ReMo::map_t::const_iterator i = remos.begin();
		remos.end() != i;
		++i)
	{
		result.insert(i->first);
	}

	return result;
}

struct RangeSorter
{
	bool operator()(ReMo::ptr_t lhs, ReMo::ptr_t rhs)
	{
		return lhs->range.start < rhs->range.start;
	}
};

bool ReMoRange::intersects(const ReMoRange & rhs) const
{
	if (start <= rhs.start)
	{
		return end > rhs.start;
	}
	else
	{
		return rhs.end > start;
	}
}

size_t
ReMoRange::get_size() const
{
	if (end < start)
	{
		throw std::logic_error( "Range end < range start" );
	}

	return end - start + 1; //end & start are inclusive 
}

double
ReMoRange::overlap(const ReMoRange & rhs) const
{
	const int overlap_start = std::max(start, rhs.start);
	const int overlap_end = std::min(end, rhs.end);
	const size_t overlap_length = std::max(0, overlap_end - overlap_start);
	return double(overlap_length) / double(get_size());
}

seq_t
ReMoBundle::get_sequence(const std::string & seq_id, bool masked)
{
	ReMo::map_t::iterator s = remos.find(seq_id);
	if (remos.end() == s)
	{
		throw BIO_MAKE_STRING("Cannot find sequence: " << seq_id);
	}

	seq_t result;
	ReMoRange last_range(std::numeric_limits<int>::min(), std::numeric_limits<int>::min());

	//the ranges should be in order but we may have to sort just in case...
	s->second.sort(RangeSorter());

	for (ReMo::list_t::const_iterator i = s->second.begin();
		s->second.end() != i;
		++i)
	{
		//check our assumptions
		BOOST_ASSERT(last_range.start <= i->get()->range.start);

		const seq_t & sequence = masked ? i->get()->masked_sequence : i->get()->unmasked_sequence;
		unsigned overlap =
			(last_range.end >= i->get()->range.start)
				? last_range.end - i->get()->range.start + 1
				: 0;
		result.append(sequence, overlap, sequence.size() - overlap);

		last_range = i->get()->range;
	}

	return result;
}


void
ReMoExtraction::serialise(const boost::filesystem::path & file) const
{
	boost::filesystem::ofstream stream(file, std::ios::binary);
    boost::archive::binary_oarchive(stream) << *this;
}

ReMoExtraction::ptr_t
ReMoExtraction::deserialise(const boost::filesystem::path & file)
{
	ptr_t result(new ReMoExtraction);

	boost::filesystem::ifstream stream(file, std::ios::binary);
	boost::archive::binary_iarchive(stream) >> *result;

	return result;
}


bool
parse_sequence_id(const std::string & value, ReMoSequenceId & id)
{
	//the gene id
	const std::string::size_type s1 = value.find(' ');
	if (std::string::npos == s1)
	{
		return false;
	}
	id.gene_id = value.substr(0, s1);

	//the transcript id
	const std::string::size_type s2 = value.find(" (", s1 + 1);
	if (std::string::npos == s2)
	{
		return false;
	}
	id.transcript_id = value.substr(s1 + 1, s2 - s1 - 1);

	//the species
	const std::string::size_type s3 = value.find(")", s2 + 2);
	if (std::string::npos == s3)
	{
		return false;
	}
	id.species = value.substr(s2 + 2, s3 - s2 - 2);

	//the index
	id.index = s3 + 1 < value.size() ? atoi(value.substr(s3 + 1).c_str()) : 0;

	return true;
}

std::string
get_gene_id(const std::string & sequence_group_name)
{
	ReMoSequenceId id;
	if (! parse_sequence_id(sequence_group_name, id))
	{
		throw BIO_MAKE_STRING("Could not parse sequence group name: \"" << sequence_group_name << "\"");
	}
	return id.gene_id;
}


void
ReMoExtraction::build_gene_remo_map(GeneReMoMap & map) const
{
	map.clear();

	for (ReMoSequenceGroup::list_t::const_iterator sg = sequence_groups.begin();
		sequence_groups.end() != sg;
		++sg)
	{
		for (ReMoBundle::map_t::const_iterator rb = sg->get()->remo_bundles.begin();
			sg->get()->remo_bundles.end() != rb;
			++rb)
		{
			ReMoBundle::id_set_t ids = rb->second->get_sequence_ids();
			for (ReMoBundle::id_set_t::const_iterator id = ids.begin();
				ids.end() != id;
				++id)
			{
				ReMoSequenceId seq_id;
				if (parse_sequence_id(*id, seq_id))
				{
					map[seq_id.gene_id].insert(rb->second);
				}
			}
		}
	}
}

std::ostream &
operator<<(std::ostream & os, ReMoLocation loc)
{
	switch(loc)
	{
	case REMO_LOC_UPSTREAM: os << "upstream"; break;
	case REMO_LOC_DOWNSTREAM: os << "downstream"; break;
	case REMO_LOC_GENEREGION: os << "gene_region"; break;
	case REMO_LOC_UNDEFINED: os << "<undefined location>"; break;
	default:
		throw std::logic_error( "Unknown ReMoLocation" );
	}
	return os;
}

std::ostream &
operator<<(std::ostream & os, ReMoSpecies species)
{
	switch(species)
	{
	case REMO_SPECIES_CHICK: os << "chick"; break;
	case REMO_SPECIES_CHIMP: os << "chimp"; break;
	case REMO_SPECIES_CIONA: os << "ciona"; break;
	case REMO_SPECIES_COW: os << "cow"; break;
	case REMO_SPECIES_DOG: os << "dog"; break;
	case REMO_SPECIES_FLY: os << "fly"; break;
	case REMO_SPECIES_HUMAN: os << "human"; break;
	case REMO_SPECIES_FUGU: os << "fugu"; break;
	case REMO_SPECIES_MOUSE: os << "mouse"; break;
	case REMO_SPECIES_OPOSSUM: os << "opossum"; break;
	case REMO_SPECIES_RAT: os << "rat"; break;
	case REMO_SPECIES_TETRAODON: os << "tetraodon"; break;
	case REMO_SPECIES_XENOPUS: os << "xenopus"; break;
	case REMO_SPECIES_ZEBRAFISH: os << "zebrafish"; break;
	case REMO_SPECIES_UNKNOWN: os << "<unknown species>"; break;
	default:
		throw std::logic_error( "Unknown ReMoSpecies" );
	}
	return os;
}


ReMoSequence::ptr_t
ReMoSequenceGroup::get_sequence_for(const std::string & id)
{
	for (ReMoSequence::list_t::const_iterator s = sequences.begin();
		sequences.end() != s;
		++s)
	{
		if (s->get()->id == id)
		{
			return *s;
		}
	}
	throw BIO_MAKE_STRING("Could not find sequence \"" << id << "\""); 
}

std::ostream &
operator<<(std::ostream & os, const ReMoRange & range)
{
	os << '[' << range.start << ',' << range.end << ']';
	return os;
}


std::ostream &
operator<<(std::ostream & os, const ReMoSubSequence & sub_sequence)
{
	os << sub_sequence.id << ',' << sub_sequence.range;
	return os;
}

std::ostream &
operator<<(std::ostream & os, const ReMoSequence & sequence)
{
	os
        << sequence.id
        << ", length=" << sequence.length
        << ", " << sequence.location;
	if (sequence.has_position)
	{
		os << ", has_position=no";
	}
	else
	{
		os << ", position=" << sequence.position;
	}

	os
        << ", " << sequence.species
        << ", exons=";

	for (ReMoExon::list_t::const_iterator e = sequence.exons.begin();
		sequence.exons.end() != e;
		++e)
	{
		os << e->range << ";";
	}

	return os;
}

std::ostream &
operator<<(std::ostream & os, const ReMo & remo)
{
	os
		<< remo.range << ','
		<< remo.conservation << ','
		<< remo.repeat_ratio << ','
		<< remo.belief;
	return os;
}

std::ostream &
operator<<(std::ostream & os, const ReMoBundle & bundle)
{
	os << "Centre: " << bundle.centre_sequence << '\n';
	for (ReMo::map_t::const_iterator r = bundle.remos.begin();
		bundle.remos.end() != r;
		++r)
	{
		os << r->first << '\n';
		for (ReMo::list_t::const_iterator r2 = r->second.begin();
			r->second.end() != r2;
			++r2)
		{
			os << **r2 << '\n';
		}
	}

	return os;
}

std::ostream &
operator<<(std::ostream & os, const ReMoSequenceGroup & group)
{
	os << "********************** Sequences **********************\n";
	for (ReMoSequence::list_t::const_iterator s = group.sequences.begin();
		group.sequences.end() != s;
		++s)
	{
		os << **s << '\n';
	}

	os << "\n********************** Bundles **********************\n";
	for (ReMoBundle::map_t::const_iterator r = group.remo_bundles.begin();
		group.remo_bundles.end() != r;
		++r)
	{
		//os << r->first << '\n';
		os << *(r->second) << '\n';
	}

	return os;
}



BIO_NS_END
