
#ifndef BIOBASE_PSSM_H_
#define BIOBASE_PSSM_H_

#include "bio/defs.h"
#include "bio/matrix.h"
#include "bio/site.h"

BIO_NS_START




struct BiobasePssmKeyTransformer
{
	typedef TableLink result_type;

	TableLink operator()( TableLink link ) const;
	TableLink operator()( Matrix * matrix ) const;
	TableLink operator()( Matrix::ptr_t matrix ) const;
	TableLink operator()( Matrix::map_t::value_type matrix ) const;
	TableLink operator()( Site * site ) const;
	TableLink operator()( Site::ptr_t site ) const;
	TableLink operator()( Site::map_t::value_type site ) const;
};




BIO_NS_END

#endif //BIOBASE_PSSM_H_


