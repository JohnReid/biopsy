/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_DEFS_H_
#define BIOPSY_DEFS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

//need to be included first
#ifndef DONT_USE_BOOST_SERIALIZATION
#ifdef _MSC_VER
# pragma warning( push )
# pragma warning( disable : 4099 )
# pragma warning( disable : 4996 )
#endif //_MSC_VER
#  include <boost/archive/text_oarchive.hpp>
#  include <boost/archive/text_iarchive.hpp>
#  include <boost/archive/binary_oarchive.hpp>
#  include <boost/archive/binary_iarchive.hpp>
#  include <boost/serialization/map.hpp>
#  include <boost/serialization/set.hpp>
#  include <boost/serialization/shared_ptr.hpp>
#  include <boost/serialization/vector.hpp>
#  include <boost/serialization/list.hpp>
#ifdef _MSC_VER
# pragma warning( pop )
#endif //_MSC_VER
#endif

#include <boost/config.hpp>

#include <boost/test/utils/wrap_stringstream.hpp>
#define BIOPSY_MAKE_STRING(x) (boost::wrap_stringstream().ref() << x).str()


#include <boost/iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/operators.hpp>
#include <boost/array.hpp>
#include <boost/regex.hpp>
#include <boost/multi_array.hpp>
#include <boost/functional/hash.hpp>

#include <string>
#include <vector>
#include <numeric>

#include <float.h>



namespace biopsy {

typedef std::string sequence;
typedef std::vector< sequence > sequence_vec;
typedef boost::shared_ptr< sequence_vec > sequence_vec_ptr;
void append_sequence( sequence_vec & v, const sequence & seq );
void append_sequences( sequence_vec & v, const sequence_vec & seqs );

typedef std::vector< std::string > string_vec;
typedef boost::shared_ptr< string_vec > string_vec_ptr;



} //namespace biopsy




#if __cplusplus > 199711L // C++11
# define BIOPSY_ISNAN( x ) std::isnan( x )
#else // not C++11
# ifdef _MSC_VER
#  define BIOPSY_ISNAN( x ) _isnan( x )
# else
#  define BIOPSY_ISNAN( x ) std::isnan( x )
# endif
#endif //__cplusplus > 199711L



#endif //BIOPSY_DEFS_H_

