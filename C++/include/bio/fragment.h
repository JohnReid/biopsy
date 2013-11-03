#ifndef BIOBASE_FRAGMENT_H_
#define BIOBASE_FRAGMENT_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/database_ref.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>

BIO_NS_START




struct Fragment : BiobaseTableEntry
{
	typedef boost::shared_ptr<Fragment> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	TableLink accession_number;
	FactorLinkList factor_links;
	seq_t sequence;
	TableLinkVec genes;

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
		ar & factor_links;
		ar & sequence;
		ar & genes;
    }		 
};





BIO_NS_END


#endif //BIOBASE_FRAGMENT_H_



