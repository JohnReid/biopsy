/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/background_models.h"
#include "bio/hmm_dna.h"
#include "bio/environment.h"

#include <fstream>

BIO_NS_START


void
gen_sequence_from_random_species_dna_hmm(seq_t & seq, size_t seq_length)
{
	DnaHmmOrderNumStateMap::singleton().gen_sequence_from_random_hmm(seq, seq_length);
}

DnaModel::ptr_t
DnaModelBroker::get_hmm_model(unsigned num_states, unsigned order)
{
	if (! DnaHmmOrderNumStateMap::singleton().contains_model(num_states, order))
	{
		throw BIO_MAKE_STRING("Do not have model with " << num_states << " states and order " << order);
	}
	return DnaHmmOrderNumStateMap::singleton().get_model_ptr( num_states, order );
}

DnaModel::ptr_t
DnaModelBroker::get_nmica_model()
{
	return DnaModel::ptr_t();
}




BIO_NS_END

