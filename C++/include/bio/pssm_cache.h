#ifndef BIO_PSSM_CACHE_H_
#define BIO_PSSM_CACHE_H_

#include "bio/defs.h"
#include "bio/cache.h"
#include "bio/singleton.h"
#include "bio/pssm.h"
#include "bio/biobase_pssm.h"
#include "bio/unary_compose.h"

BIO_NS_START



struct PssmMaker
	: std::unary_function< TableLink, boost::shared_ptr< Pssm > >
{
	boost::shared_ptr< Pssm > operator()( const TableLink & link ) const;
};


struct PssmCache
	: unary_compose<
		Dereference< Pssm >,
		unary_compose<
			Cache< PssmMaker >,
			BiobasePssmKeyTransformer
		>
	>
	, Singleton< PssmCache >
{
};





BIO_NS_END

#endif //BIO_PSSM_CACHE_H_
