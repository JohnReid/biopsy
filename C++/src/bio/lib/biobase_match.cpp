/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/biobase_match.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"


BIO_NS_START

Pssm make_pssm( TableLink link )
{
	switch( link.table_id )
	{
	case MATRIX_DATA:
		return make_pssm( BiobaseDb::singleton().get_entry< MATRIX_DATA >( link ) );
	case SITE_DATA:
		return make_pssm( BiobaseDb::singleton().get_entry< SITE_DATA >( link ) );
	default:
		throw std::logic_error( BIO_MAKE_STRING( "Cannot make_pssm() for data of type: " << link.table_id ) );
	}
}

Pssm make_pssm( const Matrix * matrix )
{
	return matrix->pssm;
}

Pssm make_pssm( Matrix::ptr_t matrix )
{
	return make_pssm( matrix.get() );
}

bool is_matchable(const Matrix * matrix)
{
	return matrix->pssm.size() > 0;
}

bool is_matchable(const Matrix::ptr_t matrix)
{
	return is_matchable(matrix.get());
}


Pssm make_pssm( const Site * site )
{
	return make_pssm_from_iupac(site->sequence.begin(), site->sequence.end());
}

Pssm make_pssm( const Site::ptr_t site )
{
	return make_pssm(site.get());
}



bool is_matchable(const Site * site)
{
	return site->sequence.size() > 0;
}

bool is_matchable(const Site::ptr_t site)
{
	return is_matchable(site.get());
}

BIO_NS_END
