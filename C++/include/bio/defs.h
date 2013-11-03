

#ifndef BIO_DEFS_H_
#define BIO_DEFS_H_

#define FUSION_MAX_VECTOR_SIZE 19

//
// Turn off uninitialised warnings in boost code.
//
#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wuninitialized"
#endif //__GNUC__

#include <boost/version.hpp>
#include <boost/config.hpp>


#ifndef BIO_NO_NAMESPACES
# define BIO_NS bio
# define BIO_NS_START namespace BIO_NS {
# define BIO_NS_END }
# define USING_BIO_NS using namespace BIO_NS;
#else //BIO_NO_NAMESPACES
# define BIO_NS
# define BIO_NS_START
# define BIO_NS_END
# define USING_BIO_NS
#endif //BIO_NO_NAMESPACES


//need to be included first
#ifndef DONT_USE_BOOST_SERIALIZATION
# ifdef _MSC_VER
#  pragma warning( push )
#  pragma warning( disable : 4099 )
#  pragma warning( disable : 4996 )
# endif //_MSC_VER
#  include <boost/archive/text_oarchive.hpp>
#  include <boost/archive/text_iarchive.hpp>
#  include <boost/archive/binary_oarchive.hpp>
#  include <boost/archive/binary_iarchive.hpp>
#  include <boost/serialization/map.hpp>
#  include <boost/serialization/set.hpp>
#  include <boost/serialization/shared_ptr.hpp>
#  include <boost/serialization/vector.hpp>
#  include <boost/serialization/list.hpp>
# ifdef _MSC_VER
#  pragma warning( pop )
# endif //_MSC_VER
#endif


#ifdef _MSC_VER
# pragma warning(push)
#  pragma warning(disable : 4312)
#  pragma warning(disable : 4311)
//    warning C4503: 'boost::mpl::vector3<T0,T1,T2>' : decorated name length exceeded, name was truncated
#  pragma warning(disable : 4503)
#endif //_MSC_VER
#include <boost/any.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/array.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/concept_check.hpp>
#include <boost/config.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/local_time/local_time_io.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/function.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/adapted.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/operators.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/random.hpp>
#include <boost/range/adaptor/indirected.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm_ext/is_sorted.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/numeric.hpp>
#include <boost/regex.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/static_assert.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/utils/wrap_stringstream.hpp>
#include <boost/timer.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility.hpp>
#include <boost/version.hpp>
#ifdef _MSC_VER
# pragma warning(pop)
#endif //_MSC_VER


//
// Turn warnings back on
//
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif //__GNUC__



#if BOOST_VERSION >= 104400
# if BOOST_VERSION >= 104900
#  define BIO_FPC_NS ::boost::test_tools ///< Floating point comparison namespace
# else
#  define BIO_FPC_NS ::boost::math::fpc ///< Floating point comparison namespace
# endif
#else
# define BIO_FPC_NS ::boost::test_tools ///< Floating point comparison namespace
#endif


//
// Conditional define to ensure works with older filesystem v2.
//
#if BOOST_FILESYSTEM_VERSION < 2
#pragma error("Can only use Boost.Filesystem version 3 or better. You will need Boost 1.44 or newer: see Boost documentation.")
#endif
#define _BOOST_FS_NATIVE string




#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define BIO_MAKE_STRING(x) (boost::wrap_stringstream().ref() << x).str()

#include <float.h>
#ifdef sun
# include <ieeefp.h>
#endif //sun
#ifdef _MSC_VER
# define BIO_ISNAN(x) _isnan(x)
# define BIO_FINITE(x) _finite(x)
# define STRICMP _stricmp
#else
#ifdef __APPLE__
# define BIO_ISNAN(x) std::isnan(x)
#else
# define BIO_ISNAN(x) isnan(x)
#endif
# define BIO_FINITE(x) finite(x)
# define STRICMP strcasecmp
#endif


#include <string>
/** Converts a string to upper case in place. */
#define BIO_TO_UPPER(s) (std::transform((s).begin(), (s).end(), (s).begin(), (int(*)(int)) toupper), s)




namespace boost { namespace program_options {
    //forward decl
    class options_description;
} }




BIO_NS_START

typedef float float_t;
typedef double prob_t;

/** Returns the argument converted to upper case. */
inline
std::string bio_to_upper(const std::string & str)
{
    std::string result(str);
    return BIO_TO_UPPER(result);
}



template<
    typename InIt,
    typename OutIt>
inline
OutIt copy_at_most(InIt begin, InIt end, unsigned n, OutIt dest)
{
    for ( ; end != begin && 0 != n; ++begin, ++dest, --n)
    {
        *dest = *begin;
    }

    return dest;
}


/**
 * Dereference a shared pointer. Used to work around some problems with boost::spirit in versions 1.52 and above
 * due to result_of change.
 */
template< typename T >
T &
deref_shared_ptr( boost::shared_ptr< T > & ptr ) {
    return *ptr;
}


BIO_NS_END

#ifdef min
# undef min
#endif //min

#ifdef max
# undef max
#endif //max



#endif //BIO_DEFS_H_

