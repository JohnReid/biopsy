/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/pssm_cache.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/biobase_match.h"


BIO_NS_START



boost::shared_ptr< Pssm >
PssmMaker::operator()( const TableLink & link ) const
{
	switch (link.table_id)
	{
	case MATRIX_DATA:
		return boost::shared_ptr< Pssm >( new Pssm( make_pssm( BiobaseDb::singleton().get_entry< MATRIX_DATA >( link ) ) ) );

	case SITE_DATA:
		return boost::shared_ptr< Pssm >( new Pssm( make_pssm( BiobaseDb::singleton().get_entry< SITE_DATA >( link ) ) ) );

	default:
		throw BIO_MAKE_STRING( link << ": Cannot make a pssm from this type of biobase entry" );
	}
}

BIO_NS_END
