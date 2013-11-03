#ifndef BIOBASE_DATABASE_REF_H_
#define BIOBASE_DATABASE_REF_H_

#include "bio/defs.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <string>
#include <vector>

BIO_NS_START

//forward decl
struct TableLink;


enum Database {
	UNKNOWN_DB,
	AFFY_PROBE_DB,
	BKL_DB,
	CELLLINE_DB,
	DIP_DB,
	EMBL_DB,
	ENSEMBL_DB,
	ENTREZ_GENE_DB,
	ENTREZ_PROTEIN_DB,
	EPD_DB,
	FLYBASE_DB,
	INPARANOID_DB,
	JASPAR_DB,
	MGI_DB,
	PATHO_DB,
	PDB_DB,
	PIR_DB,
	PROSITE_DB,
	REFSEQ_DB,
	RGD_DB,
	RSNP_DB,
	SGD_DB,
	SMART_DB,
	SWISSPROT_DB,
	TAIR_DB,
	TRANSCOMPEL_DB,
	TRANSFAC_DB,
	TRANSPATH_DB,
	TRANSPRO_DB,
	UNIGENE_DB,
	WORMBASE_DB,
	ZFIN_DB,
};
Database parse_database_string(const std::string & database);
const char * get_database_name( Database database );
std::ostream &
operator<<( std::ostream & os, Database database );


struct db_ref
{
	typedef boost::shared_ptr< db_ref > ptr;

	Database db; /**< The database the reference refers to. */
	std::string table; /**< The table in the database. */
	int acc; /**< The accession number in the table. */

	db_ref();

	db_ref( 
		Database db,
		const std::string & table,
		int acc );

	int compare( const db_ref & rhs ) const;

	bool operator==( const db_ref & rhs ) const;
	bool operator<( const db_ref & rhs ) const;

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & db;
        ar & table;
        ar & acc;
    }		 
};
typedef db_ref DatabaseRef;
typedef std::vector< DatabaseRef > DatabaseRefVec;


/** Parses acc as database reference. Returns UNKNOWN_DB if can't. */
db_ref
try_to_parse_db_ref( const std::string & acc );

/** Converts a transfac table link to a db_ref. */
db_ref
db_ref_from_transfac_table_link( const TableLink & link );

/** Gets the url for the reference. */
std::string url_for( const db_ref & ref );

/** Converts a db_ref to a transfac table link if possible. */
TableLink
transfac_table_link_from_db_ref( const db_ref & ref );

/** Parses acc as database reference. Throws exception if can't. */
db_ref 
parse_db_ref( const std::string & acc );

db_ref 
parse_db_ref_as( const std::string & acc, Database db );

/** Parses database references as held in transfac. */
std::pair< bool, db_ref >
parse_transfac_db_ref( const std::string & ref_string );

/** Parses database reference entire lines as held in transfac. */
std::pair< bool, db_ref >
parse_transfac_db_ref_line( const std::string & line );


std::size_t hash_value( const db_ref & ref );

std::ostream &
operator<<( std::ostream & os, const db_ref & ref );



BIO_NS_END

#endif //BIOBASE_DATABASE_REF_H_



