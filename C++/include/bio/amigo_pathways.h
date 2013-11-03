#ifndef BIO_AMIGO_PATHWAYS_H_
#define BIO_AMIGO_PATHWAYS_H_


/**
 * Some data structures for pathway data.
 */

#include "bio/defs.h"
#include "bio/database_ref.h"

#include <string>
#include <map>
#include <vector>



BIO_NS_START


extern const char * pathway_xml_files [];


struct AmigoPathway {
	std::string name;
	DatabaseRefVec database_refs;
	bool contains_database_ref(const DatabaseRef & database_ref) const;
	std::string get_name() const { return name; }
};
typedef boost::shared_ptr<AmigoPathway> AmigoPathwayPtr;
typedef std::vector<AmigoPathwayPtr> AmigoPathwayVec;
typedef std::map<std::string, AmigoPathwayPtr> AmigoPathwayMap;

AmigoPathwayPtr parse_pathway_xml( const std::string & xml_filename );

struct AmigoPathwaySet {
	template <class FilenameIt>
	void
	parse_xml_files(FilenameIt begin, FilenameIt end)
	{
		for ( ; begin != end; ++begin) {
			AmigoPathwayPtr pathway = parse_pathway_xml(*begin);
			pathways[pathway->name] = pathway;
		}
	}

	template <class InsIt>
	void find_pathways_containing(const DatabaseRef & database_ref, InsIt insert_it) const
	{
		for(AmigoPathwayMap::const_iterator i = pathways.begin(); i != pathways.end(); ++i) {
			if (i->second->contains_database_ref(database_ref)) {
				*insert_it++ = i->second;
			}
		}
	}

	void parse_default_xml_files();

	AmigoPathwayMap pathways;
};

extern AmigoPathwaySet pathways;

BIO_NS_END


#endif //BIO_AMIGO_PATHWAYS_H_



