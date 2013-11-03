#ifndef BIOBASE_MOLECULE_H_
#define BIOBASE_MOLECULE_H_

#include "bio/defs.h"
#include "bio/common.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>

BIO_NS_START




struct Molecule : BiobaseTableEntry
{
	typedef boost::shared_ptr<Molecule> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	virtual ~Molecule() { }

	TableLink accession_number;
	TableLinkVec pathways;
	TableLinkVec super_families; //molecular super families
	DatabaseRefVec database_refs;

	TableLink get_link() const { return accession_number; }
	std::string get_name() const { return accession_number.get_text(); }
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
        ar & pathways;
        ar & super_families;
        ar & database_refs;
    }		 
};





BIO_NS_END


#endif //BIOBASE_MOLECULE_H_



