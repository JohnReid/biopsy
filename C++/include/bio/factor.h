#ifndef BIOBASE_FACTOR_H_
#define BIOBASE_FACTOR_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/database_ref.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>
#include <set>

BIO_NS_START




struct Factor : BiobaseTableEntry
{
	typedef boost::shared_ptr<Factor> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;
	typedef std::set< std::string > synonym_set_t;
	typedef std::vector< std::string > taxonomies_t;
	
	enum type {
        family_type,
        family_mod_type,
        isogroup_type,
        isogroup_mod_type,
        basic_type,
        basic_mod_type,
        complex_type,
        complex_mod_type,
        miRNA_type,
        pre_miRNA_type,
		unknown_type
	};

	Factor() : _type( unknown_type ) { }
	virtual ~Factor() { }

    TableLink accession_number;
	DatabaseRefVec database_refs;
	std::string name;
	synonym_set_t synonyms;
	taxonomies_t taxonomies;
	TableLinkVec matrices;
	TableLink gene;
	type _type;
	TableLinkVec subunits;
	TableLinkVec complexes;
	TableLinkVec sub_families;
	TableLinkVec super_families;

	static bool parse_db_type( Database type );

	TableLink get_link() const { return accession_number; }
	std::string get_name() const { return name; }
	std::string get_description() const { return get_name(); }
	bool is_synonym(std::string n) const
	{
		BIO_TO_UPPER(n);

		//is the factor's name the same as the factor we're looking for?
		if (bio_to_upper(name) == n)
		{
			return true;
		}

		//do any synonyms match?
		for (synonym_set_t::const_iterator s = synonyms.begin();
			synonyms.end() != s;
			++s)
		{
			if (bio_to_upper(*s) == n)
			{
				return true;
			}
		}

		return false;
	}

	bool is_vertebrate() const
	{
		return std::find(taxonomies.begin(), taxonomies.end(), "vertebrata") != taxonomies.end();
	}

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & accession_number;
        ar & name;
        ar & database_refs;
		ar & synonyms;
		ar & taxonomies;
		ar & matrices;
		ar & gene;
		ar & _type;
		ar & subunits;
		ar & complexes;
		ar & sub_families;
		ar & super_families;
    }		 
};

Factor::type parse_factor_type( const std::string & s );


struct BiobaseFactorKeyTransformer
{
	typedef TableLink result_type;

	const TableLink & operator()( const TableLink & link ) const;
	const TableLink & operator()( const Factor * ptr ) const;
	const TableLink & operator()( const Factor::ptr_t ptr ) const;
	const TableLink & operator()( const Factor::map_t::value_type factor ) const;
	const TableLink & operator()( const FactorLink & factor_link ) const;
	const TableLink & operator()( const FactorLinkPtr & factor_link ) const;
};

struct CaseInsensitiveCmp
{
	bool operator()(const TableLink & s1, const TableLink & s2) const;
};

/** Used when building SVG. */
struct FactorInfo
{
	FactorInfo() : score(0.0) { }
	float_t score;
	std::set< size_t > hits;
};

typedef std::map<TableLink, FactorInfo, CaseInsensitiveCmp> factor_scores_map_t;


BIO_NS_END


#endif //BIOBASE_SITE_H_



