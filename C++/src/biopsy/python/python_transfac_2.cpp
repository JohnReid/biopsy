/**
@file

Copyright John Reid 2006

*/

#include <boost/python.hpp>
#include "biopsy/transfac.h"
#include "biopsy/db_ref.h"
#include "biopsy/python.h"

#include <bio/biobase_filter.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
#include <bio/factor.h>
#include <bio/gene.h>
#include <bio/matrix.h>
#include <bio/database_ref.h>


using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;


namespace biopsy {


namespace detail {

using namespace BIO_NS;

template< typename Container >
std::string
comma_separated_string( Container const & c )
{
	ostringstream os;
	os << "[";
	BOOST_FOREACH( typename Container::value_type const & v, c ) os << v << ",";
	os << "]";
	return os.str();
}

typedef boost::shared_ptr< TableLink > table_link_ptr;
inline table_link_ptr make_table_link( const std::string & s ) {
	return table_link_ptr(
		new TableLink( BIO_NS::parse_table_link_accession_number( s ) ) );
}

bool is_link_known( const TableLink & link ) { return UNKNOWN_DATA != link.table_id && link.entry_idx >= 0; }

object get_table_link_entry( const TableLink & link );

table_link_ptr
table_link_from_db_ref( const db_ref & ref )
{
	const TableLink link = transfac_table_link_from_db_ref( ref );
	return UNKNOWN_DATA == link.table_id ? table_link_ptr() : table_link_ptr( new TableLink( link ) );
}

struct db_ref_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs(db_ref const& r)
	{
		return boost::python::make_tuple();
	}

	static
    boost::python::tuple
    getstate(const db_ref& r)
    {
		return boost::python::make_tuple(int(r.db), r.table, r.acc);
    }

    static
    void
    setstate(db_ref& r, boost::python::tuple state)
    {
		r.db = BIO_NS::Database(int(extract<int>(state[0])));
		r.table = extract<std::string>(state[1]);
		r.acc = extract<int>(state[2]);
    }
};

struct table_link_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs(TableLink const& r)
	{
		return boost::python::make_tuple();
	}

	static
    boost::python::tuple
    getstate(const TableLink& r)
    {
		return boost::python::make_tuple(int(r.table_id), r.entry_idx);
    }

    static
    void
    setstate(TableLink& r, boost::python::tuple state)
    {
		r.table_id = BIO_NS::TransData(int(extract<int>(state[0])));
		r.entry_idx = extract<int>(state[1]);
    }
};

std::size_t table_link_hash_value( const TableLink & link )
{
	std::size_t seed = 0;
	boost::hash_combine( seed, int( link.table_id ) );
	boost::hash_combine( seed, link.entry_idx );

	return seed;
}


} // namespace detail


//
// export
//
void
export_transfac_2()
{

	//std::cout << detail::get_transfac_entry_from_acc< BIO_NS::FACTOR_DATA >( 28 ) << std::endl;

	using boost::python::arg;

	class_< BIO_NS::BiobasePssmFilter >(
		"PssmFilter",
		boost::python::init<
			bool,
			const std::string &,
			const std::string &
		>(
			(
				arg( "use_consensus_sequences" ) = true,
				arg( "species_filter" ) = "V",
				arg( "name_regex_pattern" ) = "."
			),
			"Filters pssms (sites and matrices) based on certain attributes"
		)
	).def(
		"all_pssms",
		&BIO_NS::BiobasePssmFilter::get_all_pssms_filter )
	.staticmethod( "all_pssms" )
	;



	//
	// TransData
	//
	using BIO_NS::TransData;
	enum_< TransData >( "trans_data" )
	.value( "cell", BIO_NS::CELL_DATA )
	.value( "matrix", BIO_NS::MATRIX_DATA )
	.value( "site", BIO_NS::SITE_DATA )
	.value( "factor", BIO_NS::FACTOR_DATA )
	.value( "fragment", BIO_NS::FRAGMENT_DATA )
	.value( "gene", BIO_NS::GENE_DATA )
	.value( "reference", BIO_NS::REFERENCE_DATA )
	.value( "compel", BIO_NS::COMPEL_DATA )
	.value( "evidence", BIO_NS::EVIDENCE_DATA )
	.value( "pathway", BIO_NS::PATHWAY_DATA )
	.value( "molecule", BIO_NS::MOLECULE_DATA )
	.value( "reaction", BIO_NS::REACTION_DATA )
	.value( "s", BIO_NS::S_DATA )
	.value( "unknown_data", BIO_NS::UNKNOWN_DATA )
	;


	//
	// TableLink
	//
	using BIO_NS::TableLink;
	class_< TableLink >(
		"TableLink",
		"A link to entry in a TRANSFAC database",
		init< TransData, int >(
			( arg( "table" ) = BIO_NS::UNKNOWN_DATA, arg( "entry" ) = 0 )
		)
	)
	.def(
		init< std::string >() )
	.def(
		"__init__",
		make_constructor( detail::table_link_from_db_ref ) )
	.def_pickle( detail::table_link_pickle_suite() )
	.def(
		"__eq__",
		&TableLink::operator== )
	.def(
		"__hash__",
		detail::table_link_hash_value )
	.add_property(
		"known",
		detail::is_link_known,
		"Is the link valid?" )
	.def_readonly(
		"db",
		&TableLink::table_id,
		"The database the link is to." )
	.def_readonly(
		"acc",
		&TableLink::entry_idx,
		"The accession number." )
	.def(
		"__repr__",
		as_string< TableLink >,
		"The link as a string." )
	.add_property(
		"entry",
		detail::get_table_link_entry,
		"The entry in TRANSFAC for the link." )
	.add_property(
		"url",
		&TableLink::get_url,
		"The url of the entry." )
	.add_property(
		"name",
		&TableLink::get_name,
		"The name of the entry." )
	.def(
		"as_db_ref",
		BIO_NS::db_ref_from_transfac_table_link,
		"The table link as a DbRef." )
	.def(
		"from_db_ref",
		BIO_NS::transfac_table_link_from_db_ref,
		"Create a table link from a DbRef." )
	.staticmethod( "from_db_ref" )
	;
	implicitly_convertible< std::string, TableLink >();
	register_ptr_to_python< detail::table_link_ptr >();
	using BIO_NS::TableLinkVec;
	class_< TableLinkVec >(
		"TableLinkVec",
		"A sequence of links to TRANSFAC entries" )
        .def( container_suite< TableLinkVec >())
		.def( "__repr__", detail::comma_separated_string< TableLinkVec >, "A string representation" );
    ;


	//
	// FactorLink
	//
	using BIO_NS::FactorLink;
	class_< FactorLink >(
		"FactorLink",
		"A link to factor in a TRANSFAC database",
		no_init )
	.def_readonly(
		"link",
		&FactorLink::link,
		"TableLink for the factor" )
	.def_readonly(
		"quality",
		&FactorLink::quality,
		"Quality of the association" )
	.def_readonly(
		"name",
		&FactorLink::name,
		"Name of the factor" )
	.def_readonly(
		"species",
		&FactorLink::species,
		"The species of the factor" )
	.def_readonly(
		"cellular_source",
		&FactorLink::cellular_source )
	.def(
		"__str__",
		&FactorLink::get_text,
		"String representation" )
	;
	register_ptr_to_python< BIO_NS::FactorLinkPtr >();
	using BIO_NS::FactorLinkList;
	class_< FactorLinkList >(
		"FactorLinkList",
		"A sequence of links to factor links" )
    .def( container_suite< FactorLinkList >())
    ;


	//
	// Identifier
	//
	using BIO_NS::Identifier;
	class_< Identifier >(
		"Identifier",
		"Identifies a factor",
		no_init )
	.def_readonly(
		"species_group",
		&Identifier::species_group,
		"Species" )
	.def_readonly(
		"factor",
		&Identifier::factor,
		"Factor" )
	.def_readonly(
		"discriminator",
		&Identifier::discriminating_extension,
		"Discriminates" )
	.def(
		"__str__",
		&Identifier::get_text,
		"String representation" )
	;


	//
	// Database
	//
	using BIO_NS::Database;
	enum_< Database >( "db" )
	.value( "affy", BIO_NS::AFFY_PROBE_DB )
	.value( "bkl", BIO_NS::BKL_DB )
	.value( "dip", BIO_NS::DIP_DB )
	.value( "embl", BIO_NS::EMBL_DB )
	.value( "ensembl", BIO_NS::ENSEMBL_DB )
	.value( "entrez_gene", BIO_NS::ENTREZ_GENE_DB )
	.value( "entrez_protein", BIO_NS::ENTREZ_PROTEIN_DB )
	.value( "epd", BIO_NS::EPD_DB )
	.value( "flybase", BIO_NS::FLYBASE_DB )
	.value( "inparanoid", BIO_NS::INPARANOID_DB )
	.value( "mgi", BIO_NS::MGI_DB )
	.value( "patho", BIO_NS::PATHO_DB )
	.value( "pdb", BIO_NS::PDB_DB )
	.value( "pir", BIO_NS::PIR_DB )
	.value( "refseq", BIO_NS::REFSEQ_DB )
	.value( "rgd", BIO_NS::RGD_DB )
	.value( "rsnp", BIO_NS::RSNP_DB )
	.value( "sgd", BIO_NS::SGD_DB )
	.value( "smart", BIO_NS::SMART_DB )
	.value( "swissprot", BIO_NS::SWISSPROT_DB )
	.value( "tair", BIO_NS::TAIR_DB )
	.value( "transcompel", BIO_NS::TRANSCOMPEL_DB )
	.value( "transfac", BIO_NS::TRANSFAC_DB )
	.value( "transpath", BIO_NS::TRANSPATH_DB )
	.value( "unigene", BIO_NS::UNIGENE_DB )
	.value( "wormbase", BIO_NS::WORMBASE_DB )
	.value( "zfin", BIO_NS::ZFIN_DB )
	.value( "unknown_db", BIO_NS::UNKNOWN_DB )
	;
	def(
		"database_name",
		BIO_NS::get_database_name,
		"The name of a database" );


	//
	// Database Ref
	//
	using BIO_NS::db_ref;
	class_<
		db_ref
	>(
		"DbRef",
		"Common format to store database references",
		init< BIO_NS::Database, std::string, int >()
	)
	.def( init<>() )
	.def_pickle( detail::db_ref_pickle_suite() )
	.def(
		"__cmp__",
		&db_ref::compare,
		"Compare 2 objects" )
	.def(
		"__hash__",
		( std::size_t (*)( const db_ref & ) ) BIO_NS::hash_value,
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
	.add_property(
		"url",
		( std::string ( * )( const db_ref & ) ) BIO_NS::url_for,
		"The url for this database reference" )
	.def(
		"parse_as",
		BIO_NS::parse_db_ref_as,
		( arg( "acc" ), arg( "database" ) ),
		"Parse a database reference as a particular database type" )
	.staticmethod( "parse_as" )
	.def(
		"try_to_parse",
		BIO_NS::try_to_parse_db_ref,
		"Try to parse a database reference. Return an UNKNOWN_DB if can't" )
	.staticmethod( "try_to_parse" )
	.def(
		"parse",
		BIO_NS::parse_db_ref,
		"Parse a database reference. Raise an exception if can't" )
	.staticmethod( "parse" )
	;
	using BIO_NS::DatabaseRefVec;
	class_< DatabaseRefVec >(
		"DbRefVec",
		"A sequence of database references" )
        .def( container_suite< DatabaseRefVec >() )
		.def( "__repr__", detail::comma_separated_string< DatabaseRefVec >, "A string representation" );
    ;

	export_transfac_3();

}




} //namespace biopsy


