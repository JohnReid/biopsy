
#ifndef BIO_MATRIX_MATCH_H_
#define BIO_MATRIX_MATCH_H_

#include "bio/defs.h"
#include "bio/common.h"	

#include <boost/filesystem/path.hpp>

#include <map>	


BIO_NS_START

struct MatrixMatch
{
	MatrixMatch();

	typedef std::map<TableLink, MatrixMatch> map_t;

	TableLink accession_number;
	float_t threshold;
};



void
parse_match_thresholds(
	const boost::filesystem::path & match_thresholds_file,
	MatrixMatch::map_t & match_map);


const MatrixMatch::map_t & get_min_fp_match_map();
const MatrixMatch::map_t & get_min_fn_match_map();
const MatrixMatch::map_t & get_min_sum_match_map();

BIO_NS_END


#endif //BIO_MATRIX_MATCH_H_

