
#ifndef BIO_TRANSCRIPTION_FACTOR_H_
#define BIO_TRANSCRIPTION_FACTOR_H_

#include "bio/defs.h"
#include "bio/binding_hit.h"

BIO_NS_START



/**
A transcription factor.
*/
struct TF
{
	typedef boost::shared_ptr< TF > ptr_t;
	typedef BindingHit< TF > hit_t;
	typedef BindingHitSet< TF >::type hit_set_t;
	typedef boost::shared_ptr< hit_set_t > hit_set_ptr_t;

	virtual ~TF();

	virtual std::string get_name() const = 0;
};







BIO_NS_END

#endif //BIO_TRANSCRIPTION_FACTOR_H_
