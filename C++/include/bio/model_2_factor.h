
#ifndef BIO_MODEL_2_FACTOR_H_
#define BIO_MODEL_2_FACTOR_H_

#include "bio/defs.h"
#include "bio/binding_model.h"
#include "bio/transcription_factor.h"

BIO_NS_START




void model_hits_2_factor_hits(
	const BindingModel::hit_set_t & model_hits,
	TF::hit_set_t & factor_hits );




BIO_NS_END

#endif //BIO_MODEL_2_FACTOR_H_

