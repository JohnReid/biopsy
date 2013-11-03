#ifndef BIOBASE_SITE_H_
#define BIOBASE_SITE_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/database_ref.h"
#include "bio/pathway_associations.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>

BIO_NS_START

struct Site : BiobaseTablePssmEntry
{
	typedef boost::shared_ptr<Site> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	TableLink accession_number;
	Identifier id;
	seq_t sequence;
	std::string description;
	FactorLinkList factor_links;
	DatabaseRefVec database_refs;
	std::string reference_point;
	int start_position;
	int end_position;

	Site() : start_position(0), end_position(0) { }
	TableLink get_link() const { return accession_number; }
	std::string get_name() const { return id.get_text(); }
	std::string get_description() const { return description; }
	size_t get_size() const { return sequence.size(); }
	const FactorLinkList & get_factors() const { return factor_links; }

	static bool parse_db_type( Database type );

	int get_transfac_gene_accession() const;

	TableLink get_most_significant_pathway() const
	{
		return PathwayAssociations::singleton().get_most_significant_pathway_for( accession_number );
	};
	
private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & id;
        ar & accession_number;
        ar & sequence;
        ar & description;
        ar & factor_links;
        ar & database_refs;
        ar & reference_point;
        ar & start_position;
        ar & end_position;
    }		 
};




BIO_NS_END


#endif //BIOBASE_SITE_H_



