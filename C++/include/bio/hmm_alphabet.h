#ifndef BIO_HMM_ALPHABET_H_
#define BIO_HMM_ALPHABET_H_

#include "bio/defs.h"



BIO_NS_START


/** Specialise this for particular alphabet's details. */
template <class Alphabet>
struct AlphabetTraits {
	/**
	static size_t get_size();
	static size_t get_index(Alphabet a);
	static Alphabet get_symbol(size_t index);
	*/
};


BIO_NS_END

#endif //BIO_HMM_ALPHABET_H_
