#ifndef BIO_PSSM_LIKELIHOOD_H_
#define BIO_PSSM_LIKELIHOOD_H_


#include "bio/defs.h"
#include "bio/pssm.h"
#include "bio/biobase_likelihoods.h"



BIO_NS_START

/** Calculates the likelihoods of certain biobase scores given sequences that match the pssm. */ 
void
pssm_likelihood(
	const Pssm & pssm,
	const size_t max_map_size,
	BiobaseLikelihoods & result,
	float_t pseudo_count = 1.0,
	bool verbose = false);

BIO_NS_END


#endif //BIO_PSSM_LIKELIHOOD_H_
