/**
@file

Copyright John Reid 2006

*/

#include "biopsy/defs.h"
#include "biopsy/transfac.h"
#include "biopsy/sequence.h"

#include <bio/matrix_match.h>
#include <bio/biobase_filter.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
USING_BIO_NS;

namespace biopsy
{

bio::BiobaseTablePssmEntry * 
get_transfac_pssm_entry( const std::string & pssm_name )
{
	if( is_transfac_pssm( pssm_name ) )
	{
		TableLink link = parse_table_link_accession_number( pssm_name );
		switch( link.table_id )
		{
		case SITE_DATA: 
		case MATRIX_DATA:
			return BiobaseDb::singleton().get_pssm_entry( link );
		default:
			break;
		}
	}
	return 0;
}


std::string
get_transfac_pssm_accession( 
	const std::string & pssm_name )
{
	typedef std::map< std::string, std::string > acc_map;

	static acc_map _acc_map;
	static bool _inited = false;

	if( ! _inited )
	{
		BiobasePssmFilter filter = BiobasePssmFilter::get_all_pssms_filter();

		BOOST_FOREACH( 
			const Matrix::map_t::value_type & pssm, 
			get_matrices( filter ) )
		{
			_acc_map[ pssm.second->get_name() ] = BIO_MAKE_STRING( pssm.second->accession_number );
		}
		
		BOOST_FOREACH( 
			const Site::map_t::value_type & pssm, 
			get_sites( filter ) )
		{
			_acc_map[ pssm.second->get_name() ] = BIO_MAKE_STRING( pssm.second->accession_number );
		}
		
		_inited = true;
	}
	acc_map::iterator i = _acc_map.find( pssm_name );
	if( _acc_map.end() == i )
	{
		throw 
			std::logic_error(
				BIO_MAKE_STRING(
					"Could not find biobase pssm with name: " << pssm_name ) );
	}
	return i->second;
}


std::string
get_transfac_pssm_name( 
	const std::string & pssm )
{
	return BiobaseDb::singleton().get_pssm_entry( parse_table_link_accession_number( pssm ) )->get_name();
}


string_vec_ptr
get_transfac_pssm_sequences(
	const std::string & pssm )
{
	USING_BIO_NS;

	string_vec_ptr result( new string_vec );
	const TableLink link = parse_table_link_accession_number( pssm );
	switch( link.table_id )
	{
	case SITE_DATA: 
		return result;

	case MATRIX_DATA:
		{
			Matrix * matrix = BiobaseDb::singleton().get_entry< MATRIX_DATA >( link );
			BOOST_FOREACH( AlignDescPtr align_desc, matrix->align_descs )
			{
				if( is_known_sequence()( align_desc->sequence ) )
				{
					result->push_back(
						align_desc->positive_orientation
							? align_desc->sequence
							: biopsy::reverse_complement( align_desc->sequence ) );
				}
			}
		}
		return result;

	default:
		throw std::invalid_argument( BIOPSY_MAKE_STRING( "Biobase entry is not a pssm: " << link ) );
	}
}


double
get_transfac_pssm_min_fp_threshold(
	const std::string & pssm )
{
	const TableLink link = parse_table_link_accession_number( pssm );
	MatrixMatch::map_t::const_iterator i = get_min_fp_match_map().find( link );
	if( get_min_fp_match_map().end() == i )
	{
		throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find min fp threshold for: " << link ) );
	}
	return i->second.threshold;
}



double
get_transfac_pssm_min_fn_threshold(
	const std::string & pssm )
{
	const TableLink link = parse_table_link_accession_number( pssm );
	MatrixMatch::map_t::const_iterator i = get_min_fn_match_map().find( link );
	if( get_min_fn_match_map().end() == i )
	{
		throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find min fn threshold for: " << link ) );
	}
	return i->second.threshold;
}

double
get_transfac_pssm_min_sum_threshold(
	const std::string & pssm )
{
	const TableLink link = parse_table_link_accession_number( pssm );
	MatrixMatch::map_t::const_iterator i = get_min_sum_match_map().find( link );
	if( get_min_sum_match_map().end() == i )
	{
		throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find min sum(fp,fn) threshold for: " << link ) );
	}
	return i->second.threshold;
}





string_vec_ptr
get_transfac_pssm_accessions( const BIO_NS::BiobasePssmFilter & filter )
{
	USING_BIO_NS;

	string_vec_ptr result( new string_vec );

	BOOST_FOREACH( const Matrix::map_t::value_type & p, get_matrices( filter ) )
	{
		result->push_back( BIOPSY_MAKE_STRING( p.first ) );
	}
	
	BOOST_FOREACH( const Site::map_t::value_type & p, get_sites( filter ) )
	{
		result->push_back( BIOPSY_MAKE_STRING( p.first ) );
	}
	
	return result;
}

bio::BiobasePssmFilter
get_default_transfac_pssm_filter()
{
	return bio::BiobasePssmFilter();
}

string_vec_ptr
get_factors_for_pssm( const std::string & pssm_acc )
{
	string_vec_ptr result( new string_vec );

	BiobaseTablePssmEntry * pssm = BiobaseDb::singleton().get_pssm_entry( parse_table_link_accession_number( pssm_acc ) );
	BOOST_FOREACH( FactorLinkPtr f, pssm->get_factors() )
	{
		result->push_back( BIOPSY_MAKE_STRING( f->link ) );
	}

	return result;
}


string_vec_ptr
get_pssms_for_factor( const std::string & factor_acc )
{
	string_vec_ptr result( new string_vec );

	const TableLink factor_link = parse_table_link_accession_number( factor_acc );

	//get all the matrices and sites that list this factor...
	bio::BiobasePssmFilter filter = bio::BiobasePssmFilter::get_all_pssms_filter();
	BOOST_FOREACH( Matrix::map_t::value_type v, get_matrices( filter ) )
	{
		BOOST_FOREACH( FactorLinkPtr f, v.second->get_factors() )
		{
			if( f->link == factor_link )
			{
				result->push_back( BIOPSY_MAKE_STRING( v.first ) );
			}
		}
	}
	BOOST_FOREACH( Site::map_t::value_type v, get_sites( filter ) )
	{
		BOOST_FOREACH( FactorLinkPtr f, v.second->get_factors() )
		{
			if( f->link == factor_link )
			{
				result->push_back( BIOPSY_MAKE_STRING( v.first ) );
			}
		}
	}

	//make the matrices and sites unique...
	result->erase( std::unique( result->begin(), result->end() ), result->end() );

	return result;
}

namespace {

template< typename BiobaseMap >
string_vec_ptr
get_accessions( const BiobaseMap & map)
{
	USING_BIO_NS;

	string_vec_ptr result( new string_vec );
	BOOST_FOREACH( typename BiobaseMap::value_type v, map )
	{
		result->push_back( BIOPSY_MAKE_STRING( v.first ) );
	}

	return result;
}

} //namespace

string_vec_ptr get_transfac_matrices( ) { return get_accessions( BiobaseDb::singleton().get_matrices() ); }
string_vec_ptr get_transfac_sites( ) { return get_accessions( BiobaseDb::singleton().get_sites() ); }
string_vec_ptr get_transfac_factors( ) { return get_accessions( BiobaseDb::singleton().get_factors() ); }
string_vec_ptr get_transfac_fragments( ) { return get_accessions( BiobaseDb::singleton().get_fragments() ); }
string_vec_ptr get_transfac_genes( ) { return get_accessions( BiobaseDb::singleton().get_genes() ); }


string_vec_ptr
get_factors_for_fragment( const std::string & fragment_acc )
{
	USING_BIO_NS;

	string_vec_ptr result( new string_vec );

	Fragment * fragment = BiobaseDb::singleton().get_entry< FRAGMENT_DATA >( parse_table_link_accession_number( fragment_acc ) );
	//get the factors for this fragment
	BOOST_FOREACH( FactorLinkPtr f, fragment->factor_links )
	{
		result->push_back( BIOPSY_MAKE_STRING( f->link ) );
	}

	return result;
}

string_vec_ptr
get_fragments_for_factor( const std::string & factor_acc )
{
	USING_BIO_NS;

	const TableLink factor = parse_table_link_accession_number( factor_acc );

	string_vec_ptr result( new string_vec );
	BOOST_FOREACH( Fragment::map_t::value_type v, BiobaseDb::singleton().get_fragments() )
	{
		BOOST_FOREACH( FactorLinkPtr f, v.second->factor_links )
		{
			if( f->link == factor )
			{
				result->push_back( BIOPSY_MAKE_STRING( v.first ) );
				break;
			}
		}
	}

	return result;
}

std::string
get_fragment_sequence( const std::string & fragment_acc )
{
	USING_BIO_NS;

	string_vec_ptr result( new string_vec );

	Fragment * fragment = BiobaseDb::singleton().get_entry< FRAGMENT_DATA >( parse_table_link_accession_number( fragment_acc ) );
	return fragment->sequence;
}



} //namespace biopsy

