#ifndef BIOBASE_PATHWAY_H_
#define BIOBASE_PATHWAY_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/database_ref.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>
#include <iostream>

BIO_NS_START



enum PathwayType
{
	CHAIN_PW,
	EVIDENCE_CHAIN_PW,
	PATHWAY_PW,
	UNKNOWN_PW,
};
inline
std::ostream &
operator<<(std::ostream & os, PathwayType pathway_type)
{
	switch (pathway_type) {
		case CHAIN_PW: os << "chain"; break;
		case EVIDENCE_CHAIN_PW: os << "evidence chain"; break;
		case PATHWAY_PW: os << "pathway"; break;
		case UNKNOWN_PW: os << "<unknown pathway type>"; break;
	}
	return os;
}


/** Represents an entry in the Biobase Pathway table. */
struct Pathway : BiobaseTableEntry {
	typedef boost::shared_ptr<Pathway> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	virtual ~Pathway() { }
	
    TableLink accession_number;
	PathwayType pathway_type;
	TableLinkVec super_families; //pathway super families
	std::string name;

	Pathway() : pathway_type(UNKNOWN_PW) { }

	TableLink get_link() const { return accession_number; }
	std::string get_name() const { return name; }
	std::string get_description() const { return get_name(); }

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & accession_number;
        ar & pathway_type;
        ar & super_families;
        ar & name;
    }		 
};

/** Pathways GK deems interesting. */
extern TableLinkVec interesting_pathways;

/** Test if pathway is interesting. */
struct IsInterestingPathway
{
	bool operator()(const TableLink & link) const;
};






BIO_NS_END


#endif //BIOBASE_PATHWAY_H_



