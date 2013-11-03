/**
@file

Copyright John Reid 2006

*/

#include "biopsy/db_ref.h"
#include "biopsy/python.h"



using namespace boost;
using namespace boost::python;
using namespace std;


namespace biopsy {

namespace detail {

BIO_NS::db_ref::ptr parse_db_ref( const std::string & acc )
{
	return BIO_NS::db_ref::ptr( new BIO_NS::db_ref( biopsy::parse_db_ref( acc ) ) );
}

} //namespace detail

void export_db_ref()
{
	class_<
		db_ref
	>(
		"DbRef",
		"Common format to store database references",
		init< BIO_NS::Database, std::string, int >()
	)
	.def(
		"__init__",
		make_constructor( detail::parse_db_ref ) )
	.def(
		"__cmp__",
		&db_ref::compare,
		"Compare 2 objects" )
	.def(
		"__hash__",
		( std::size_t (*)( const db_ref & ) ) hash_value,
		"Hash value for use in dict e.g." )
	.def(
		"__repr__",
		as_string< db_ref >,
		"String representation" )
	.def_readonly(
		"db",
		&db_ref::db,
		"The database the reference is in" )
	.def_readonly(
		"table",
		&db_ref::table,
		"The table in the database" )
	.def_readonly(
		"acc",
		&db_ref::acc,
		"The accession number of the database" )
	;
	class_< DatabaseRefVec >(
		"DbRefVec",
		"A sequence of database references" )
        .def( container_suite< DatabaseRefVec >() )
    ;

}



} //namespace biopsy
