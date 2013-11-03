#ifndef BIOBASE_COMMON_H_
#define BIOBASE_COMMON_H_


/**
 * Some common data structures for biobase data.
 */

#include "bio/defs.h"
#include "bio/database_ref.h"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/access.hpp>
#include <boost/operators.hpp>

#include <list>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

BIO_NS_START


enum TransData
{
	//Transfac
	CELL_DATA,
	MATRIX_DATA,
	SITE_DATA,
	FACTOR_DATA,
	FRAGMENT_DATA,
	GENE_DATA,
	REFERENCE_DATA,

	//Transcompel
	COMPEL_DATA,
	EVIDENCE_DATA,

	//Transpath
	PATHWAY_DATA,
	MOLECULE_DATA,
	REACTION_DATA,

	//unknown
	S_DATA,
	UNKNOWN_DATA,
};
TransData parse_biobase_table( const std::string & table_name );
const char * biobase_table_to_string( TransData table );
Database database_from_trans_data_type( TransData table );

inline
std::ostream &
operator<<(std::ostream & os, TransData table)
{
	switch (table)
	{
		case CELL_DATA: os << "Cell"; break;
		case COMPEL_DATA: os << "Compel"; break;
		case EVIDENCE_DATA: os << "Evidence"; break;
		case FACTOR_DATA: os << "Factor"; break;
		case FRAGMENT_DATA: os << "Fragment"; break;
		case GENE_DATA: os << "Gene"; break;
		case MATRIX_DATA: os << "Matrix"; break;
		case MOLECULE_DATA: os << "Molecule"; break;
		case PATHWAY_DATA: os << "Pathway"; break;
		case REACTION_DATA: os << "Reaction"; break;
		case REFERENCE_DATA: os << "Reference"; break;
		case S_DATA: os << "CompelSite"; break;
		case SITE_DATA: os << "Site"; break;
		case UNKNOWN_DATA: os << "<unknown>"; break;
		default:
			throw std::logic_error( "Unknown biobase table" );
	}
	return os;
}


/** A link to an entry in some table. Like "R00130". */
struct TableLink
	: boost::equality_comparable<TableLink>
	, boost::less_than_comparable<TableLink>
{
	TableLink() : table_id(UNKNOWN_DATA), entry_idx(0) { }
	TableLink(
		TransData table_id,
		int entry_idx);
	explicit TableLink( const std::string & s );

	TransData table_id;
	int entry_idx;

	std::string get_text() const;
	std::string get_url() const;
	bool operator==(const TableLink & rhs) const;
	bool operator<(const TableLink & rhs) const;
	std::string get_name() const;

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & table_id;
        ar & entry_idx;
    }		 
};
typedef std::vector<TableLink> TableLinkVec;
TableLink parse_table_link_accession_number(const std::string & accession_number);

std::ostream &
operator<<(std::ostream & os, const TableLink & tl);

struct FactorLink
	: boost::equality_comparable< FactorLink >
	, boost::less_than_comparable< FactorLink >
{
	TableLink link;
	int quality;
	std::string name;
	std::vector<std::string> species;
	std::string cellular_source;
	bool sites_included;

	std::string get_text() const;

	bool operator==( const FactorLink & rhs ) const;
	bool operator<( const FactorLink & rhs ) const;

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & link;
        ar & quality;
        ar & name;
        ar & species;
        ar & cellular_source;
        ar & sites_included;
    }		 
};
typedef boost::shared_ptr< FactorLink > FactorLinkPtr;
typedef std::vector< FactorLinkPtr > FactorLinkList;
bool contains_factor( const FactorLinkList & factor_links, TableLink factor ); 

inline
std::ostream &
operator<<(std::ostream & os, const FactorLink & fl) {
	os << fl.link << " " << fl.name << " ";
	BOOST_FOREACH( const std::string & species, fl.species )
		{os << species << " ";}
	return os;
}

inline
std::ostream &
operator<<(std::ostream & os, const FactorLinkList & fl) {
	for (FactorLinkList::const_iterator i = fl.begin(); i != fl.end(); ++i) {
		os << **i << "; ";
	}
	return os;
}



enum Species
{
	MOUSE_SPECIES,
	RAT_SPECIES,
	HUMAN_SPECIES,
	DOG_SPECIES,
	XENOPUS_SPECIES,
	TETRAODON_SPECIES,
	COW_SPECIES,
	CHICKEN_SPECIES,
	FUGU_SPECIES,
	UNKNOWN_SPECIES
};

inline
std::ostream &
operator<<(std::ostream & os, Species species)
{
	switch(species) {
		case MOUSE_SPECIES: os << "mouse"; break;
		case RAT_SPECIES: os << "rat"; break;
		case HUMAN_SPECIES: os << "human"; break;
		case DOG_SPECIES: os << "dog"; break;
		case XENOPUS_SPECIES: os << "xenopus"; break;
		case TETRAODON_SPECIES: os << "tetraodon"; break;
		case COW_SPECIES: os << "cow"; break;
		case CHICKEN_SPECIES: os << "chicken"; break;
		case FUGU_SPECIES: os << "fugu"; break;
		case UNKNOWN_SPECIES: os << "<unknown species>"; break;
		default: throw std::logic_error( "Undefined species" );
	}
	return os;
}


struct Identifier {
	std::string species_group;
	std::string factor;
	std::string discriminating_extension;

	const std::string get_text() const;

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & species_group;
        ar & factor;
        ar & discriminating_extension;
    }		 
};

inline
std::ostream &
operator<<(std::ostream & os, const Identifier & id) {
	return os << id.species_group << "$" << id.factor << "_" << id.discriminating_extension;
}

/** An interface with some common methods for entries in Biobase tables. */
struct BiobaseTableEntry
{
	virtual ~BiobaseTableEntry() { }

	virtual TableLink get_link() const = 0;
	virtual std::string get_name() const = 0;
	virtual std::string get_description() const = 0;

	bool operator==( const BiobaseTableEntry & rhs ) const { return get_link() == rhs.get_link(); }
	bool operator<( const BiobaseTableEntry & rhs ) const { return get_link() < rhs.get_link(); }

};

/** An interface with some common methods for entries in Biobase tables matrix and site. */
struct BiobaseTablePssmEntry
: BiobaseTableEntry
{
	virtual ~BiobaseTablePssmEntry() { }

	virtual size_t get_size() const = 0;
	virtual const FactorLinkList & get_factors() const = 0;
	bool is_vertebrate() const;
};



BIO_NS_END


#endif //BIOBASE_COMMON_H_



