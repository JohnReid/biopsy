#ifndef BIO_EVIDENCE_H_
#define BIO_EVIDENCE_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/database_ref.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>

BIO_NS_START




struct Evidence : BiobaseTableEntry
{
	typedef boost::shared_ptr<Evidence> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	virtual ~Evidence() { }

	TableLink accession_number;
    TableLink composite_element; //the composite element the evidence is for
	DatabaseRefVec database_refs;

	static bool parse_db_type( Database type );

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
        ar & composite_element;
		ar & database_refs;
    }		 
};





BIO_NS_END


#endif //BIO_EVIDENCE_H_



