#ifndef BIOBASE_NORMALISER_TRANSFAC_H_
#define BIOBASE_NORMALISER_TRANSFAC_H_

#include "bio/matrix.h"
#include "bio/site.h"
#include "bio/ott_match_normaliser.h"

BIO_NS_START


template <class MatrixIt, class SeqIt>
size_t
cache_normalisations(
	MatrixIt matrix_begin,
	MatrixIt matrix_end,
	SeqIt seq_begin,
	SeqIt seq_end,
	BIO_NS::OttNormaliserCache & cache,
	size_t num_quanta)
{
	size_t num_normalised = 0;
	//for each matrix
	for ( ; matrix_begin != matrix_end; ++matrix_begin) {
		const std::string name = matrix_begin->second->id.get_text();
		OttNormaliser::ptr_t ott_normaliser(new OttNormaliser(make_pssm(matrix_begin->second), seq_begin, seq_end, num_quanta));
		cache.put_normaliser(name, ott_normaliser);
		++num_normalised;
	}
	return num_normalised;
}



BIO_NS_END


#endif //BIOBASE_NORMALISER_TRANSFAC_H_
