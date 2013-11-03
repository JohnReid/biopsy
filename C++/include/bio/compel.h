#ifndef BIO_COMPEL_H_
#define BIO_COMPEL_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <map>

BIO_NS_START


namespace detail {
template< typename T, typename Str >
T
safe_lexical_cast( const Str & s )
{
	try
	{
		return boost::lexical_cast< T >( s );
	}
	catch( const std::exception & e )
	{
		throw std::logic_error( BIO_MAKE_STRING( "\"" << s << "\": " << e.what() ) );
	}
}
} //namespace detail


/**
Represents the data in a line such as:
BS  -64 to -55; NF-kappaB or p50:RelA-p65; s00087.
in the compel.dat TransCompel table.
*/
struct CompelBindingSite {
	typedef std::vector< CompelBindingSite > vec;

	int start;
	int end;
	std::string factor; /**< E.g. "NF-kappaB or p50:RelA-p65" */
	TableLink site_link; /**< E.g. s00087 */

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & start;
        ar & end;
        ar & factor;
        ar & site_link;
    }		 
};


enum CompelType {
	COMPEL_SYNERGISM,		//synergism
	COMPEL_ANTAGONISM,		//antagonism
	COMPEL_TYPE_UNKNOWN,	//not known
};


/**
An entry in the TransCompel compel table.
*/
struct Compel : BiobaseTableEntry
{
	typedef boost::shared_ptr<Compel> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	virtual ~Compel() { }

	TableLink accession_number;
    Identifier id;
	TableLink gene;
	seq_t sequence;
	int begin;
	int end;
	CompelBindingSite::vec binding_sites;
	CompelType type;
	DatabaseRefVec database_refs;
	std::string comment;
	TableLinkVec evidences;

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
        ar & id;
		ar & gene;
		ar & sequence;
		ar & begin;
		ar & end;
		ar & binding_sites;
		ar & type;
		ar & database_refs;
		ar & comment;
        ar & evidences;
    }		 
};





BIO_NS_END


#endif //BIO_COMPEL_H_



