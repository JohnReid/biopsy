#ifndef BIO_REMO_H_
#define BIO_REMO_H_

#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/common.h"

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/serialization/access.hpp>

#include <antlr/TokenStreamSelector.hpp>

#include <map>
#include <set>



BIO_NS_START

struct Remo {
	typedef std::map<Species, seq_t> remo_map_t;
	remo_map_t map;
};
typedef boost::shared_ptr<Remo> RemoPtr;

template <class OStr>
OStr &
operator<<(OStr & os, const Remo & remo)
{
	for (Remo::remo_map_t::const_iterator i = remo.map.begin(); remo.map.end() != i; ++i) {
		os << i->first << ":" << i->second.size() << ":" << i->second << std::endl;
	}
	return os;
}

typedef std::map<std::string, RemoPtr> RemoMap;
template <class OStr>
OStr &
operator<<(OStr & os, const RemoMap & remo_map)
{
	for (RemoMap::const_iterator i = remo_map.begin(); remo_map.end() != i; ++i)
	{
		os
			<< i->first << std::endl
			<< *(i->second) << std::endl
			<< std::endl;
	}
	return os;
}


/** Some test remos. */
extern RemoMap test_remos;

/** Populate the test remos. */
void build_test_remos();

enum ReMoLocation
{
	REMO_LOC_UPSTREAM,
	REMO_LOC_DOWNSTREAM,
	REMO_LOC_GENEREGION,
	REMO_LOC_UNDEFINED,
};
std::ostream &
operator<<(std::ostream & os, ReMoLocation loc);



enum ReMoSpecies
{
	REMO_SPECIES_CHICK,
	REMO_SPECIES_CHIMP,
	REMO_SPECIES_CIONA,
	REMO_SPECIES_COW,
	REMO_SPECIES_DOG,
	REMO_SPECIES_FLY,
	REMO_SPECIES_FUGU,
	REMO_SPECIES_HUMAN,
	REMO_SPECIES_MOUSE,
	REMO_SPECIES_OPOSSUM,
	REMO_SPECIES_RAT,
	REMO_SPECIES_TETRAODON,
	REMO_SPECIES_XENOPUS,
	REMO_SPECIES_ZEBRAFISH,
	REMO_SPECIES_UNKNOWN,
};
ReMoSpecies parse_remo_species(const std::string & species_string);
std::ostream &
operator<<(std::ostream & os, ReMoSpecies loc);



struct EnsemblId
{
	EnsemblId(const std::string & id = "") : value(id) { }

	std::string value;
	unsigned number;

	bool operator==(const EnsemblId & rhs) const { return value == rhs.value; }

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & value;
		ar & number;
    }
};




/**
A range of a sequence.
*/
struct ReMoRange
{
	int start; /**< Inclusive. */
	int end; /**< Inclusive. */

	ReMoRange(int start = -1, int end = -1)
		: start(start), end(end)
	{ }

	bool operator<(const ReMoRange & rhs) const
	{
		return start < rhs.start || (start == rhs.start && end < rhs.end);
	}

	bool operator==(const ReMoRange & rhs) const
	{
		return start == rhs.start && end == rhs.end;
	}

	bool intersects(const ReMoRange & rhs) const;

	double overlap(const ReMoRange & rhs) const;

	size_t get_size() const;

	int get_centre() const
	{
		return (start + end) / 2;
	}

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & start;
        ar & end;
    }
};
std::ostream &
operator<<(std::ostream & os, const ReMoRange & range);


struct ReMoSequenceId
{
	std::string gene_id;
	std::string transcript_id;
	std::string species;
	unsigned index;
};
bool parse_sequence_id(const std::string & value, ReMoSequenceId & id);
std::string get_gene_id(const std::string & sequence_group_name);


/**
A sub-sequence in a ReMoSequence is identified by its sequence and its position in that sequence.
*/
struct ReMoSubSequence
{
	std::string id;
	ReMoRange range;

	ReMoSubSequence(const std::string & id = "", const ReMoRange & range = ReMoRange())
		: id(id), range(range)
	{ }

	bool operator<(const ReMoSubSequence & rhs) const
	{
		return id < rhs.id || (id == rhs.id && range < rhs.range);
	}

	bool operator==(const ReMoSubSequence & rhs) const
	{
		return id == rhs.id && range == rhs.range;
	}

	typedef std::set<ReMoSubSequence> set_t;

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & id;
        ar & range;
    }
};
std::ostream &
operator<<(std::ostream & os, const ReMoSubSequence & sub_sequence);


struct ReMoExon
{
	EnsemblId id;
	ReMoRange range;

	ReMoExon(const EnsemblId & id = EnsemblId(), const ReMoRange & range = ReMoRange()) : id(id), range(range) { }

	typedef std::list<ReMoExon> list_t;

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & id;
        ar & range;
    }
};

/** A ReMoSequence defines a sequence for one species that holds remos. */
struct ReMoSequence
{
	std::string id;
	int length;
	ReMoLocation location;
	bool has_position;
	int position; /**< Position relative to TSS. */
	ReMoSpecies species;
	ReMoExon::list_t exons;

	ReMoSequence()
		: id("")
		, length(-1)
		, has_position(false)
		, position(0)
		, species(REMO_SPECIES_UNKNOWN)
	{ }

	typedef boost::shared_ptr<ReMoSequence> ptr_t;
	typedef std::list<ptr_t> list_t;

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & id;
        ar & length;
        ar & location;
        ar & has_position;
        ar & position;
        ar & species;
        ar & exons;
    }
};
std::ostream &
operator<<(std::ostream & os, const ReMoSequence & sequence);







struct ReMo
{
	ReMoRange range;
	ReMoRange target_range;
	unsigned conservation;
	unsigned repeat_ratio;
	double belief;
	seq_t masked_sequence;
	seq_t unmasked_sequence;

	ReMo()
		: conservation(0)
		, repeat_ratio(0)
		, belief(0.0)
	{ }

	const seq_t & get_sequence( bool masked = true ) const
	{
		return
			masked
				? masked_sequence
				: unmasked_sequence;
	}

	typedef boost::shared_ptr<ReMo> ptr_t;
	typedef std::list<ptr_t> list_t;
	typedef std::map<std::string, list_t> map_t; /**< Maps sequence ids to a list of remos. */

	/** Inserts the sequence for the remo (list) as one chunk into the insert iterator. */
	template < typename InsIt >
	static
	void
	copy_sequence(
		const list_t & remo_list,
		InsIt insert_it,
		bool masked = true)
	{
		for (list_t::const_iterator r = remo_list.begin();
			remo_list.end() != r;
			++r)
		{
			const seq_t & seq = r->get()->get_sequence( masked );
			std::copy( seq.begin(), seq.end(), insert_it );
		}
	}

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & range;
        ar & target_range;
        ar & conservation;
        ar & repeat_ratio;
        ar & belief;
        ar & masked_sequence;
        ar & unmasked_sequence;
    }
};
std::ostream &
operator<<(std::ostream & os, const ReMo & remo);


struct ReMoBundle
{
	ReMo::map_t remos;
	std::string centre_sequence;

	typedef boost::shared_ptr<ReMoBundle> ptr_t;
	typedef std::list<ptr_t> list_t;
	typedef std::set<ptr_t> set_t;
	typedef std::map<ReMoSubSequence, ptr_t> map_t;
	typedef std::set<std::string> id_set_t;

	ReMoBundle()
	{ }

	/** Get the id of the centre sequence of the comparison. */
	ReMoSubSequence get_centre_sequence_id() const;

	/** Get the set of sequence ids that make up this remo. */
	id_set_t get_sequence_ids() const;

	/** Get the sequence for this id. */
	seq_t get_sequence(const std::string & id, bool masked);

	const ReMo & get_centre_remo() const
	{
		ReMo::map_t::const_iterator r = remos.find(centre_sequence);
		if (remos.end() == r)
		{
			throw BIO_MAKE_STRING("Could not find centre remo: " << centre_sequence);
		}
		if (r->second.empty())
		{
			throw BIO_MAKE_STRING("Centre remo is empty: " << centre_sequence);
		}
		return *(r->second.begin()->get());
	}

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & remos;
		ar & centre_sequence;
    }
};
std::ostream &
operator<<(std::ostream & os, const ReMoBundle & bundle);



struct ReMoSequenceGroup
{
	ReMoSequence::list_t sequences;
	ReMoBundle::map_t remo_bundles;

	ReMoSequenceGroup()
	{ }

	bool operator==(const ReMoSequenceGroup & rhs) const;

	//gets the sequence for the given remo
	ReMoSequence::ptr_t get_sequence_for(const std::string & id);

	/** Get the centre sequence. */
	const ReMoSequence & get_centre_sequence() const
	{
		if (sequences.empty())
		{
			throw std::logic_error( "No sequences in group!" );
		}

		return *(sequences.begin()->get());
	};

	typedef boost::shared_ptr<ReMoSequenceGroup> ptr_t;
	typedef std::list<ptr_t> list_t;

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & sequences;
        ar & remo_bundles;
    }
};
std::ostream &
operator<<(std::ostream & os, const ReMoSequenceGroup & group);

struct ReMoExtraction
{
	ReMoSequenceGroup::list_t sequence_groups;

	typedef boost::shared_ptr<ReMoExtraction> ptr_t;
	typedef std::map<std::string, ReMoBundle::set_t> GeneReMoMap;

	void build_gene_remo_map(GeneReMoMap & map) const;

	bool operator==(const ReMoExtraction & rhs) const;

	void serialise(const boost::filesystem::path & file) const;
	static ptr_t deserialise(const boost::filesystem::path & file);

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & sequence_groups;
    }
};


ReMoExtraction::ptr_t parse_remo_extraction(boost::filesystem::path file);



struct ReMoSharedParserState
{
	antlr::TokenStreamSelector * selector;

	void push_lexer(const char * name) {
		selector->push(name);
	}
	void pop_lexer() {
		selector->pop();
	}
	void log(const char * text) {
		std::cout << text << std::endl;
	}	
};

BIO_NS_END

#endif //BIO_REMO_H_

