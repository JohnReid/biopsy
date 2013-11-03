/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_CONTAINERS_H_
#define BIOPSY_CONTAINERS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include <boost/multi_array.hpp>

#include <vector>





namespace biopsy {

typedef std::vector< double > double_vector; /**< A vector of double. */
typedef std::vector< double_vector > double_vector_vec; /**< A vector of vector of floats. */
typedef boost::multi_array< double, 2 > double_array;
typedef std::vector< double_vector_vec > double_vector_vec_vec; /**< A vector of vector of vector of doubles. */

typedef std::vector< unsigned > unsigned_vector; /**< A vector of unsigned. */
typedef std::vector< unsigned_vector > unsigned_vector_vec; /**< A vector of vector of unsigned. */
typedef boost::multi_array< unsigned, 3 > unsigned_array_3d;


} //namespace biopsy

#endif //BIOPSY_CONTAINERS_H_
