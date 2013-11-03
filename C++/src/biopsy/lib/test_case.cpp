/**
@file

Copyright John Reid 2006

*/

#include "biopsy/test_case.h"
#include "biopsy/transfac.h"

#include "bio/site_test_case.h"
#include "bio/biobase_db.h"
USING_BIO_NS



namespace biopsy {
namespace detail {

string_vec_ptr
get_pssms_for_factors( string_vec_ptr factors )
{
	string_vec_ptr result( new string_vec );

#if 0
	//get all the pssms for the factors
	BOOST_FOREACH( const std::string & f, *factors )
	{
		string_vec_ptr pssms_for_factor = get_pssms_for_factor( f );

		//add to positives
		std::copy( pssms_for_factor->begin(), pssms_for_factor->end(), std::back_inserter( *( result ) ) );
	}

	//make positives unique
	result->erase( std::unique( result->begin(), result->end() ), result->end() );
#else //implementation with a set...
	std::set< std::string > pssms;
	//get all the pssms for the factors
	BOOST_FOREACH( const std::string & f, *factors )
	{
		string_vec_ptr pssms_for_factor = get_pssms_for_factor( f );

		//add to positives
		std::copy( pssms_for_factor->begin(), pssms_for_factor->end(), std::inserter( pssms, pssms.begin() ) );
	}
	std::copy( pssms.begin(), pssms.end(), std::back_inserter( *( result ) ) );
#endif
	return result;
}


} //namespace detail

test_case::test_case( const std::string & name )
: _name( name )
, _positives( new string_vec )
, _sequences( new sequence_vec )
{
}

test_case::vec_ptr
get_transfac_site_test_cases()
{
	USING_BIO_NS;

	test_case::vec_ptr result( new test_case::vec );
	BOOST_FOREACH( SiteTestCase::ptr_t bio_tc, SiteTestCases::singleton() )
	{
		test_case::ptr biopsy_tc( new test_case( bio_tc->name ) );

		//copy the sequences into the new test case
		biopsy_tc->_sequences->push_back( bio_tc->centre_sequence );
		std::copy( 
			bio_tc->sequences.begin(), 
			bio_tc->sequences.end(), 
			std::back_inserter( *( biopsy_tc->_sequences ) ) );

		//get all the factors for this site
		string_vec_ptr factors_for_pssm = get_factors_for_pssm( BIOPSY_MAKE_STRING( bio_tc->site ) );

		//get all the pssms for these factors
		biopsy_tc->_positives = detail::get_pssms_for_factors( factors_for_pssm );

		result->push_back( biopsy_tc );

		//if( result->size() > 10 ) break;
	}

	return result;
}

test_case::vec_ptr
get_transfac_fragment_test_cases()
{
	USING_BIO_NS;

	test_case::vec_ptr result( new test_case::vec );

	BOOST_FOREACH( Fragment::map_t::value_type v, BiobaseDb::singleton().get_fragments() )
	{
		test_case::ptr biopsy_tc( new test_case( BIOPSY_MAKE_STRING( v.first ) ) );

		const std::string fragment_acc = BIOPSY_MAKE_STRING( v.first );
		string_vec_ptr factors = get_factors_for_fragment( fragment_acc );
		biopsy_tc->_positives = detail::get_pssms_for_factors( factors );

		biopsy_tc->_sequences->push_back( v.second->sequence );

		result->push_back( biopsy_tc );

		//if( result->size() > 10 ) break;
	}

	return result;
}

} //namespace biopsy
