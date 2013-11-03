/**
@file

Copyright John Reid 2006, 2013

*/

#include "biopsy/alias_lookup.h"
#include <stdexcept>
#include <boost/python/extract.hpp>

using namespace boost;
using namespace boost::assign;
using namespace std;

#define PY_STR_AS_STD( x )( std::string( boost::python::extract< const char * >( x ) ) )
namespace biopsy {

enum db_id
{
    db_unknown,
    db_entrez,
    db_ensembl,
};

/** The accession for a particular entry in a database, e.g. P3476 in UniProt. */
typedef std::string db_accession;

/** A collection of accessions. */
typedef std::vector< db_accession > accession_vec;

/** A map from one accession (of one db) to a list of accessions (in another db). */
typedef std::map< db_accession, accession_vec > accession_map;

/** A collection of maps from one particular databases' accession ids to others (indexed by their db_id). */
typedef std::map< db_id, accession_map > db_map;

/** A collection of maps for all databases (indexed by their db_id). */
typedef std::map< db_id, db_map > database_aliases;


database_aliases _aliases;

db_id get_id_for( const std::string & db_name )
{
    if( boost::iequals( db_name, "entrez" ) )
    {
        return db_entrez;
    }
    
    if( boost::iequals( db_name, "ensembl" ) 
        ||
        boost::iequals( db_name, "ens" ) )
    {
        return db_ensembl;
    }

    return db_unknown;
};

boost::python::list
lookup( 
    boost::python::str from_db_name, 
    boost::python::str from_accession, 
    boost::python::str to_db_name )
{
    std::cout << "Throwing exception" << std::endl;
    throw std::logic_error( "Error" );

    const db_id from_db = get_id_for( PY_STR_AS_STD( from_db_name ) );
    if( db_unknown == from_db )
    {
        throw std::logic_error( "Unknown 'from' database" );
    }

    const db_id to_db = get_id_for( PY_STR_AS_STD( to_db_name ) );
    if( db_unknown == to_db )
    {
        throw std::logic_error( "Unknown 'to' database" );
    }

    if( to_db == from_db )
    {
        throw std::logic_error( "'from' and 'to' databases are the same" );
        //throw std::logic_error( "'from' and 'to' databases are the same" );
    }

    const db_accession from_acc = PY_STR_AS_STD( from_accession );

    const accession_vec & to_accessions = _aliases[ from_db ][ to_db ][ from_acc ];

    boost::python::list result;
    BOOST_FOREACH( const db_accession & to_acc, to_accessions )
    {
        result.append( boost::python::str( to_acc ) );
    }

    return result;
}


void add(
    boost::python::str from_db_name,
    boost::python::str from_accession,
    boost::python::str to_db_name,
    boost::python::str to_accession )
{
    const db_id from_db = get_id_for( PY_STR_AS_STD( from_db_name ) );
    if( db_unknown == from_db )
    {
        throw std::logic_error( "Unknown 'from' database" );
    }

    const db_id to_db = get_id_for( PY_STR_AS_STD( to_db_name ) );
    if( db_unknown == to_db )
    {
        throw std::logic_error( "Unknown 'to' database" );
    }

    if( to_db == from_db )
    {
        throw std::logic_error( "'from' and 'to' databases are the same" );
    }

    const db_accession from_acc = PY_STR_AS_STD( from_accession );
    const db_accession to_acc = PY_STR_AS_STD( to_accession );

    accession_vec & to_accessions = _aliases[ from_db ][ to_db ][ from_acc ];

    if( to_accessions.end() == std::find( to_accessions.begin(), to_accessions.end(), to_acc ) )
    {
        to_accessions.push_back( to_acc );
    }
}





} //namespace biopsy


