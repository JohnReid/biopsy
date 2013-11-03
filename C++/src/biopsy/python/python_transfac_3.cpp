/**
@file

Copyright John Reid 2006

*/

#include "biopsy/python.h"
#include "biopsy/transfac.h"
#include "biopsy/db_ref.h"
#include "biopsy/sequence.h"

#include <bio/biobase_filter.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
#include <bio/factor.h>
#include <bio/gene.h>
#include <bio/matrix.h>
#include <bio/database_ref.h>

#include <iostream>

namespace boost { namespace python { namespace indexing {
template<>
struct value_traits< BIO_NS::CompelBindingSite > : value_traits< int >
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
template<>
struct value_traits< BIO_NS::PssmEntry > : value_traits< int >
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
template<>
struct value_traits< BIO_NS::AlignDescPtr > : value_traits< int >
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
} } }


using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;


#if 0
struct X {
	typedef boost::shared_ptr< X > ptr;
	std::string _s;
	X( const std::string & s = "test" ) : _s( s ) { }
};
X::ptr make_x( ) { return X::ptr( new X ); }
X::ptr make_x_from_string( const std::string & s ) { return X::ptr( new X( s ) ); }
void export_x()
{
	using namespace boost::python;
	class_< X > x( "X", no_init );
	x.def( "__init__", make_constructor( make_x ) );
	x.def( "__init__", make_constructor( make_x_from_string ) );
	x.def_readonly( "string", &X::_s );
	register_ptr_to_python< X::ptr >();
}
#endif

namespace biopsy {


namespace detail {

using namespace BIO_NS;

std::size_t table_link_hash_value( const TableLink & link );

typedef boost::shared_ptr< TableLink > table_link_ptr;

template< TransData data_type >
struct table_helper
{
	typedef DataTraits< data_type > traits;
	typedef typename traits::entry_t entry;
	typedef typename entry::ptr_t ptr;
	typedef typename entry::map_t map;
	typedef class_< entry, ptr > bp_class;

	typedef std::vector< ptr > vec;
	typedef boost::shared_ptr< vec > vec_ptr;

	static ptr get_entry( TableLink id )
	{
		typename DataTraits< data_type >::entry_t::ptr_t result = BiobaseDb::singleton().get_entry_ptr< data_type >( id );
		if( ! result )
		{
			throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find: " << id ) );
		}
		return result;
	}

	static ptr get_entry_from_acc( int acc )
	{
		return get_entry( TableLink( data_type, acc ) );
	}

	static const map & get_map();

	static vec_ptr all()
	{
		vec_ptr result( new vec );
		BOOST_FOREACH( typename map::value_type p, get_map() )
		{
			result->push_back( p.second );
		}
		return result;
	}

	static std::size_t hash_value( const entry & entry ) { return table_link_hash_value( entry.get_link() ); }
	static std::string name( const entry & entry ) { return entry.get_name(); }
	static std::string description( const entry & entry ) { return entry.get_description(); }
	static std::string link_as_string( const entry & entry ) { return BIOPSY_MAKE_STRING( entry.get_link() ); }
	static bool equals( const entry & entry_1, const entry & entry_2 ) { return entry_1 == entry_2; }


	static vec_ptr filter( const BiobasePssmFilter & filter );

	static bp_class export_class()
	{
		std::string capitalised_name = traits::get_name();
		capitalised_name[0] += ('A' - 'a');

		bp_class _c(
			capitalised_name.c_str(),
			BIOPSY_MAKE_STRING( "A " << traits::get_name() << " in TRANSFAC" ).c_str(),
			no_init );
		_c.def(
			"__init__",
			make_constructor( get_entry ) );
		_c.def(
			"__init__",
			make_constructor( get_entry_from_acc ) );
		_c.def(
			"__hash__",
			hash_value );
		_c.def(
			"__eq__",
			equals );
		_c.def(
			"__str__",
			name,
			"String representation" );
		_c.add_property(
			"name",
			name,
			"Name" );
		_c.add_property(
			"description",
			description,
			"Description" );
		_c.def(
			"__repr__",
			link_as_string,
			"String representation" );
		_c.def_readonly(
			"acc",
			&entry::accession_number,
			"Accession number in TRANSFAC" );
		_c.def(
			"all",
			all,
			BIOPSY_MAKE_STRING( "All " << traits::get_name() << " entries in TRANSFAC" ).c_str() );
		_c.staticmethod( "all" );
		class_< vec, vec_ptr >(
			BIOPSY_MAKE_STRING( capitalised_name << "Seq" ).c_str(),
			BIOPSY_MAKE_STRING( "A " << traits::get_name() << " sequence" ).c_str() )
		.def( container_suite< vec >())
		;

		return _c;
	}
};

template< > const table_helper< SITE_DATA >::map & table_helper< SITE_DATA >::get_map() {
	return BiobaseDb::singleton().get_sites();
}

template< > const table_helper< MATRIX_DATA >::map & table_helper< MATRIX_DATA >::get_map() {
	return BiobaseDb::singleton().get_matrices();
}

template< > const table_helper< FACTOR_DATA >::map & table_helper< FACTOR_DATA >::get_map() {
	return BiobaseDb::singleton().get_factors();
}

template< > const table_helper< FRAGMENT_DATA >::map & table_helper< FRAGMENT_DATA >::get_map() {
	return BiobaseDb::singleton().get_fragments();
}

template< > const table_helper< GENE_DATA >::map & table_helper< GENE_DATA >::get_map() {
	return BiobaseDb::singleton().get_genes();
}

template< > const table_helper< COMPEL_DATA >::map & table_helper< COMPEL_DATA >::get_map() {
	return BiobaseDb::singleton().get_compels();
}

template< > const table_helper< EVIDENCE_DATA >::map & table_helper< EVIDENCE_DATA >::get_map() {
	return BiobaseDb::singleton().get_evidences();
}

template< > const table_helper< PATHWAY_DATA >::map & table_helper< PATHWAY_DATA >::get_map() {
	return BiobaseDb::singleton().get_pathways();
}

template< > const table_helper< MOLECULE_DATA >::map & table_helper< MOLECULE_DATA >::get_map() {
	return BiobaseDb::singleton().get_molecules();
}

template< >
table_helper< SITE_DATA >::vec_ptr
table_helper< SITE_DATA >::filter( const BiobasePssmFilter & filter )
{
	vec_ptr result( new vec );
	BOOST_FOREACH( map::value_type p, get_sites( filter ) )
	{
		result->push_back( p.second );
	}
	return result;
}

template< >
table_helper< MATRIX_DATA >::vec_ptr
table_helper< MATRIX_DATA >::filter( const BiobasePssmFilter & filter )
{
	vec_ptr result( new vec );
	BOOST_FOREACH( map::value_type p, get_matrices( filter ) )
	{
		result->push_back( p.second );
	}
	return result;
}

Gene::ptr_t site_gene( const Site & site )
{
	return table_helper< GENE_DATA >::get_entry_from_acc( site.get_transfac_gene_accession() );
}

table_link_ptr factor_gene( const Factor & factor )
{
	return
		GENE_DATA == factor.gene.table_id
			? table_link_ptr( new TableLink( factor.gene ) )
			: table_link_ptr();
}

string_vec_ptr factor_synonyms( const Factor & factor )
{
	string_vec_ptr result( new string_vec );
	std::copy( factor.synonyms.begin(), factor.synonyms.end(), back_inserter(*result) );
	return result;
}

object get_table_link_entry( const TableLink & link )
{
	switch( link.table_id )
	{
	case MATRIX_DATA: return object( table_helper< MATRIX_DATA >::get_entry( link ) );
	case SITE_DATA: return object( table_helper< SITE_DATA >::get_entry( link ) );
	case FACTOR_DATA: return object( table_helper< FACTOR_DATA >::get_entry( link ) );
	case FRAGMENT_DATA: return object( table_helper< FRAGMENT_DATA >::get_entry( link ) );
	case GENE_DATA: return object( table_helper< GENE_DATA >::get_entry( link ) );
	case COMPEL_DATA: return object( table_helper< COMPEL_DATA >::get_entry( link ) );
	case EVIDENCE_DATA: return object( table_helper< EVIDENCE_DATA >::get_entry( link ) );
	case PATHWAY_DATA: return object( table_helper< PATHWAY_DATA >::get_entry( link ) );
	case MOLECULE_DATA: return object( table_helper< MOLECULE_DATA >::get_entry( link ) );
	default: break;
	}
	throw std::logic_error( BIOPSY_MAKE_STRING( "Cannot return entries of this type: " << link.table_id ) );
}


object
get_matrix_sequences( const Matrix & matrix )
{
	string_vec_ptr result( new string_vec );
	BOOST_FOREACH( AlignDescPtr align_desc, matrix.align_descs )
	{
		if( is_known_sequence()( align_desc->sequence ) && align_desc->sequence.size() == matrix.pssm.size() )
		{
			result->push_back(
				align_desc->positive_orientation
					? align_desc->sequence
					: biopsy::reverse_complement( align_desc->sequence ) );
		}
	}
	return object( result );
}

std::string consensus_as_string( BIO_NS::ConsensusMatrix const& cons ) {
	std::string result;
	std::copy( cons.begin(), cons.end(), std::back_inserter( result ) );
	return result;
}

boost::python::str matrix_consensus_as_string( BIO_NS::Matrix const& m ) {
	return boost::python::str( consensus_as_string( m.consensus_matrix ).c_str() );
}

std::string pssm_entry_as_string( BIO_NS::PssmEntry const& entry ) {
	USING_BIO_NS;
	return BIOPSY_MAKE_STRING( entry );
}

std::string pssm_as_string( BIO_NS::Pssm const& pssm ) {
	USING_BIO_NS;
	return BIOPSY_MAKE_STRING( pssm );
}

std::string align_desc_as_string( BIO_NS::AlignDesc const& align_desc ) {
	return BIOPSY_MAKE_STRING(
		align_desc.sequence
		<< " "
		<< align_desc.site
		<< " "
		<< (align_desc.positive_orientation ? "+ve" : "-ve")
		<< " "
		<< align_desc.start<<":"<<align_desc.length
		<< (
			BIO_NS::UNKNOWN_DATA != align_desc.secondary_site.table_id
				? BIOPSY_MAKE_STRING( " ("<<align_desc.secondary_site<<")" )
				: std::string()
		)
	);
}

char base_from_string( boost::python::str s ) {
	std::string _s = boost::python::extract< std::string >( s );
	if( _s.size() != 1 ) {
		throw std::logic_error( "String should have length 1 for conversion to base." );
	}
	return _s[0];
}

BIO_NS::float_t pssm_entry_get_count( BIO_NS::PssmEntry const& pe, boost::python::str s ) { return pe.get_count( base_from_string( s ) ); }
BIO_NS::float_t pssm_entry_get_score( BIO_NS::PssmEntry const& pe, boost::python::str s ) { return pe.get_score( base_from_string( s ) ); }
BIO_NS::float_t pssm_entry_get_freq( BIO_NS::PssmEntry const& pe, boost::python::str s, BIO_NS::float_t pseudo_count = 0.0 ) { return pe.get_freq( base_from_string( s ), pseudo_count ); }


} // namespace detail


//
// export
//
void
export_transfac_3()
{
	using boost::python::arg;

	//
	// AlignDesc
	//
	class_< BIO_NS::AlignDesc, BIO_NS::AlignDescPtr >(
		"AlignDesc",
		"How one site a Matrix was built from is aligned." )
	.def_readonly( "sequence", &BIO_NS::AlignDesc::sequence, "The sequence of the site." )
	.def_readonly( "site", &BIO_NS::AlignDesc::site, "A link to the site." )
	.def_readonly( "secondary_site", &BIO_NS::AlignDesc::secondary_site, "A link to the secondary site." )
	.def_readonly( "start", &BIO_NS::AlignDesc::start, "Where the matrix starts." )
	.def_readonly( "length", &BIO_NS::AlignDesc::length, "The length of the matrix." )
	.def_readonly( "gaps", &BIO_NS::AlignDesc::gaps, "Where there are gaps in the sequence that makes the alignment." )
	.def_readonly( "positive_orientation", &BIO_NS::AlignDesc::positive_orientation, "Which orientation the site is in." )
	.def( "__str__", detail::align_desc_as_string, "A string representation." )
	;

	//
	// int vector
	//
	class_< std::vector< int > >(
		"IntVector",
		"A sequence of integers." )
	.def( container_suite< std::vector< int > >())
    ;

	//
	// AlignDescList
	//
	class_< BIO_NS::AlignDescList >(
		"AlignDescList",
		"A sequence of AlignDesc's." )
	.def( container_suite< BIO_NS::AlignDescList >())
    ;

	//
	// PssmEntry
	//
	class_< BIO_NS::PssmEntry >(
		"PssmCountEntry",
		"The counts of particular nucleotides at one position in a PSSM."
	)
	.add_property( "num_observations", &BIO_NS::PssmEntry::get_num_observations, "The number of observations in this entry." )
	.add_property( "max_score", &BIO_NS::PssmEntry::get_max, "The highest score in this entry." )
	.add_property( "min_score", &BIO_NS::PssmEntry::get_min, "The lowest score in this entry." )
	.add_property( "conservation_information", &BIO_NS::PssmEntry::get_conservation_information, "The conservation information in this entry." )
	.def( "__str__", detail::pssm_entry_as_string, "A string representation." )
	.def( "freq", detail::pssm_entry_get_freq, ( arg( "base" ), arg( "pseudo_counts" ) = 0 ), "The frequency of the given base in this entry (assuming given pseudo-counts)." )
	.def( "count", detail::pssm_entry_get_count, "The number of observations of the given base in this entry." )
	.def( "score", detail::pssm_entry_get_score, "The score of the given base in this entry." )
    ;


	//
	// Pssm
	//
	class_< std::vector< BIO_NS::PssmEntry > >(
		"PssmEntryVec",
		"A sequence of PssmEntry's."
	)
	.def( container_suite< std::vector< BIO_NS::PssmEntry > >())
	;
	class_<
		BIO_NS::Pssm
		, boost::python::bases< std::vector< BIO_NS::PssmEntry > >
	>(
		"PssmCounts",
		"A PSSM represented as a sequence of PssmEntry's."
	)
	.def( "__str__", detail::pssm_as_string, "A string representation" )
    ;


	//
	// Matrix
	//
	using BIO_NS::Matrix;
	detail::table_helper< BIO_NS::MATRIX_DATA >::export_class()
	.def_readonly(
		"site_alignments",
		&Matrix::align_descs,
		"Descriptions of the site sequence alignments used to build this matrix." )
	.def_readonly(
		"factors",
		&Matrix::factor_links,
		"The factors for this matrix." )
	.def_readonly(
		"factor_name",
		&Matrix::factor_name,
		"The name of the factor this matrix represents." )
	.def_readonly(
		"pssm_counts",
		&Matrix::pssm,
		"The counts of observations that were used to construct this Matrix." )
	.add_property(
		"consensus",
		detail::matrix_consensus_as_string,
		"The consensus string of this Matrix." )
	.add_property(
		"size",
		&Matrix::get_size,
		"The size of this matrix's pssm." )
	.add_property(
		"pathway", 
		&Matrix::get_most_significant_pathway,
		"The most significant pathway associated with the matrix." )
	.add_property(
		"sequences", 
		detail::get_matrix_sequences,
		"The sequences that were used to construct this Matrix." )
	.def(
		"pssms",
		detail::table_helper< BIO_NS::MATRIX_DATA >::filter,
		( arg( "filter" ) = BIO_NS::BiobasePssmFilter::get_all_pssms_filter() ),
		"Get all those matrices that pass the filter." )
	.staticmethod( "pssms" )
	;




	//
	// Site
	//
	using BIO_NS::Site;
	detail::table_helper< BIO_NS::SITE_DATA >::export_class()
	.def_readonly(
		"id",
		&Site::id,
		"Identifies this site." )
	.def_readonly(
		"db_refs",
		&Site::database_refs,
		"The database references for this site." )
	.def_readonly(
		"sequence",
		&Site::sequence,
		"The sequence for this site." )
	.def_readonly(
		"reference_point",
		&Site::reference_point,
		"Where this site refers to." )
	.def_readonly(
		"start",
		&Site::start_position,
		"Start of this site." )
	.def_readonly(
		"end",
		&Site::end_position,
		"End of this site." )
	.add_property(
		"gene",
		detail::site_gene,
		"The TRANSFAC gene for this site." )
	.def_readonly(
		"factors",
		&Site::factor_links,
		"The factors for this site." )
	.add_property(
		"size",
		&Site::get_size,
		"The size of this site's sequence." )
	.add_property(
		"pathway", 
		&Site::get_most_significant_pathway,
		"The most significant pathway associated with the matrix." )
	.def(
		"pssms",
		detail::table_helper< BIO_NS::SITE_DATA >::filter,
		( arg( "filter" ) = BIO_NS::BiobasePssmFilter::get_all_pssms_filter() ),
		"Get all those sites that pass the filter." )
	.staticmethod( "pssms" )
	;



	//
	// Factor
	//
	using BIO_NS::Factor;
	enum_< Factor::type >( "factor_type" )
		.value( "family", Factor::family_type )
		.value( "isogroup", Factor::isogroup_type )
		.value( "basic", Factor::basic_type )
		.value( "complex", Factor::complex_type )
		.value( "miRNA", Factor::miRNA_type )
		.value( "unknown", Factor::unknown_type )
	;
	detail::table_helper< BIO_NS::FACTOR_DATA >::export_class()
	.def_readonly(
		"type",
		&Factor::_type,
		"The type of this factor (family, isogroup, basic, complex, miRNA or unknown)." )
	.add_property(
		"gene",
		detail::factor_gene,
		"A link to the gene the factor is the product of." )
	.add_property(
		"synonyms",
		detail::factor_synonyms,
		"Synonyms for this factor." )
	.def_readonly(
		"db_refs",
		&Factor::database_refs,
		"The database references for this factor." )
	.def_readonly(
		"taxonomies",
		&Factor::taxonomies,
		"Systematic biological classification of the species." )
	.def_readonly(
		"matrices",
		&Factor::matrices,
		"Matrix table entries providing DNA-binding profiles of the factor." )
	.def_readonly(
		"subunits",
		&Factor::subunits,
		"The subunits for this factor (if it is a complex)." )
	.def_readonly(
		"complexes",
		&Factor::complexes,
		"A list of complexes which contain this factor." )
	.def_readonly(
		"sub_families",
		&Factor::sub_families,
		"Lists entries, e.g. splice variants or family members of this isogroup/family entry." )
	.def_readonly(
		"super_families",
		&Factor::super_families,
		"Lists generic entries (isogroup or family) to which this factor belongs." )
	;


	//
	// Gene
	//
	using BIO_NS::Gene;
	detail::table_helper< BIO_NS::GENE_DATA >::export_class()
	.def_readonly(
		"db_refs",
		&Gene::database_refs,
		"The database references for this gene" )
	.def_readonly(
		"species",
		&Gene::species,
		"The gene's species" )
	;




	//
	// Fragment
	//
	using BIO_NS::Fragment;
	detail::table_helper< BIO_NS::FRAGMENT_DATA >::export_class()
	.def_readonly(
		"factors",
		&Fragment::factor_links,
		"The factors for this fragment" )
	.def_readonly(
		"sequence",
		&Fragment::sequence,
		"The sequence of this fragment" )
	.def_readonly(
		"genes",
		&Fragment::genes,
		"The genes this fragment is associated with" )
	;




	//
	// TransData
	//
	using BIO_NS::CompelType;
	enum_< CompelType >( "compel_type" )
	.value( "synergism", BIO_NS::COMPEL_SYNERGISM )
	.value( "antagonism", BIO_NS::COMPEL_ANTAGONISM )
	.value( "compel_type_unknown", BIO_NS::COMPEL_TYPE_UNKNOWN )
	;




	//
	// CompelBindingSite
	//
	using BIO_NS::CompelBindingSite;
	class_< CompelBindingSite >(
		"CompelBindingSite",
		"A binding site  in the Compel database",
		no_init
	)
	.def_readonly(
		"start",
		&CompelBindingSite::start,
		"The starting point for this compel binding site" )
	.def_readonly(
		"end",
		&CompelBindingSite::end,
		"The end point for this compel binding site" )
	.def_readonly(
		"factors",
		&CompelBindingSite::factor,
		"The factor(s) for this compel binding site" )
	.def_readonly(
		"site_link",
		&CompelBindingSite::site_link,
		"The site reference of this binding site" )
	;
	class_< CompelBindingSite::vec >(
		"CompelBindingSiteVec",
		"A sequence of binding sites for a Compel entry" )
        .def( container_suite< CompelBindingSite::vec >())
    ;






	//
	// Compel
	//
	using BIO_NS::Compel;
	detail::table_helper< BIO_NS::COMPEL_DATA >::export_class()
	.def_readonly(
		"id",
		&Compel::id,
		"Identifier of this compel entry" )
	.def_readonly(
		"gene",
		&Compel::gene,
		"The gene for this compel entry" )
	.def_readonly(
		"begin",
		&Compel::begin,
		"Identifies the beginning of the composite element within the promoter." )
	.def_readonly(
		"end",
		&Compel::end,
		"Identifies the end of the composite element within the promoter." )
	.def_readonly(
		"sequence",
		&Compel::sequence,
		"The sequence for this compel entry" )
	.def_readonly(
		"binding_sites",
		&Compel::binding_sites,
		"The binding sites for this compel entry" )
	.def_readonly(
		"type",
		&Compel::type,
		"The type of this compel entry" )
	.def_readonly(
		"database_refs",
		&Compel::database_refs,
		"The external database references for this compel entry" )
	.def_readonly(
		"comment",
		&Compel::comment,
		"The comment for this compel entry" )
	.def_readonly(
		"evidences",
		&Compel::evidences,
		"Evidences for this compel entry" )
	;




	//
	// Evidence
	//
	using BIO_NS::Evidence;
	detail::table_helper< BIO_NS::EVIDENCE_DATA >::export_class()
	.def_readonly(
		"composite",
		&Evidence::composite_element,
		"Composite element this evidence is for" )
	.def_readonly(
		"db_refs",
		&Evidence::database_refs,
		"Database references for this evidence" )
	;




	//
	// Pathway
	//
	using BIO_NS::PathwayType;
	enum_< PathwayType >( "pathway_type" )
	.value( "chain", BIO_NS::CHAIN_PW )
	.value( "evidence_chain", BIO_NS::EVIDENCE_CHAIN_PW )
	.value( "pathway", BIO_NS::PATHWAY_PW )
	.value( "unknown_pathway", BIO_NS::UNKNOWN_PW )
	;
	using BIO_NS::Pathway;
	detail::table_helper< BIO_NS::PATHWAY_DATA >::export_class()
	.def_readonly(
		"type",
		&Pathway::pathway_type,
		"Type of this pathway" )
	.def_readonly(
		"super_families",
		&Pathway::super_families,
		"Superfamilies of this pathway" )
	;




	//
	// Molecule
	//
	using BIO_NS::Molecule;
	detail::table_helper< BIO_NS::MOLECULE_DATA >::export_class()
	.def_readonly(
		"pathways",
		&Molecule::pathways,
		"Pathways of this molecule" )
	.def_readonly(
		"super_families",
		&Molecule::super_families,
		"Superfamilies of this molecule" )
	.def_readonly(
		"db_refs",
		&Molecule::database_refs,
		"The database references for this molecule" )
	;



}




} //namespace biopsy


