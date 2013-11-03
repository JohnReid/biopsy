#ifndef BIOBASE_GENE_H_
#define BIOBASE_GENE_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/database_ref.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>

BIO_NS_START




struct Gene : BiobaseTableEntry
{
	typedef boost::shared_ptr<Gene> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	virtual ~Gene() { }
	TableLink accession_number;
	DatabaseRefVec database_refs;
	std::string name;
	std::string species;

	TableLink get_link() const { return accession_number; }
	std::string get_name() const { return name; }
	std::string get_description() const { return get_name(); }
	std::string get_species() const { return species; }

	static bool parse_db_type( Database type );

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & accession_number;
		ar & database_refs;
		ar & name;
		ar & species;
    }		 
};





BIO_NS_END


#endif //BIOBASE_GENE_H_



