
#ifndef BIO_BIOBASE_TRANSCRIPTION_FACTOR_H_
#define BIO_BIOBASE_TRANSCRIPTION_FACTOR_H_

#include "bio/defs.h"
#include "bio/unary_compose.h"
#include "bio/transcription_factor.h"
#include "bio/cache.h"
#include "bio/singleton.h"
#include "bio/equivalent_factors.h"

BIO_NS_START




/**
Creates factors from biobase tablelink factor references.
*/
struct BiobaseFactor2TF
	: std::unary_function< EquivalentFactors::partition_ptr_t, TF::ptr_t >
{
	result_type operator()( argument_type factor_partition ) const;
};





/**
Caches biobase transcription factor equivalence partitions.
*/
struct BiobaseTFCache
	: unary_compose<
		Dereference< TF >,
		unary_compose<
			Cache< BiobaseFactor2TF >,
			EquivalentFactorKeyTransformer
		>
	>
	, Singleton< BiobaseTFCache >
{
};






BIO_NS_END

#endif //BIO_BIOBASE_TRANSCRIPTION_FACTOR_H_
