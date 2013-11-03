/**
Copyright John Reid 2010

@file Exposes BiFa C++ implementation to python

*/

#include <boost/python.hpp>
#include <myrrh/python/boost_range.h>
#include <myrrh/python/multi_array_to_numpy.h>
#include <boost/multi_array.hpp>
#include <boost/lambda/bind.hpp>

#include "biopsy/bifa.h"
#include "biopsy/pssm.h"
#include "biopsy/binding_hits.h"



namespace {

int major_version() { return 1; }
int minor_version() { return 0; }
std::string version() { return BIOPSY_MAKE_STRING( major_version() << "." << minor_version() ); }

using namespace biopsy::bifa;
using namespace boost::python;
using namespace myrrh::python;

/// C++ type for our PSSMs passed from python
typedef boost::multi_array_ref< double, 2 > pssm_t;



/// Can we extract type from the first element of the python sequence?
template< typename T >
bool
can_extract_type( object seq ) {
    if( len( seq ) > 0 ) {
        return extract< T >( seq[ 0 ] ).check();
    }
    return true;
}


#define BIFA_MAKE_ITERATOR_FROM_INT_SEQ( seq ) \
    make_extract_iterator< int >( py_seq_iterator( seq, 0 ) )
#define BIFA_MAKE_ITERATOR_FROM_CHAR_SEQ( seq ) \
    boost::make_transform_iterator( \
        make_extract_iterator< char >( py_seq_iterator( seq, 0 ) ), convert_char_base_to_int() )
#define BIFA_MAKE_RANGE_FROM_INT_SEQ( seq ) \
    make_boost_range< int >( seq )


/// Python interface to score_word
double
py_score_word(
    pssm_t pssm,
    object word
) {
    if( size_t( len( word ) ) < boost::size( pssm ) ) {
        throw std::logic_error( BIOPSY_MAKE_STRING( "Word is too short to be scored by PSSM: " << len( word ) << " < " << boost::size( pssm ) ) );
    }
    if( can_extract_type< int >( word ) ) {
        return score_word< DnaAlphabet >( pssm, BIFA_MAKE_ITERATOR_FROM_INT_SEQ( word ) );
    } else if( can_extract_type< char >( word ) ) {
        return score_word< DnaAlphabet >( pssm, BIFA_MAKE_ITERATOR_FROM_CHAR_SEQ( word ) );
    } else throw std::logic_error( "Cannot interpret as sequence." );
}


/// Python interface to score_one_strand
list
py_score_one_strand(
    pssm_t pssm,
    object sequence
) {
    list result;

    using namespace boost::lambda;
    using boost::ref;

    if( can_extract_type< int >( sequence ) ) {
#ifdef _WIN32
        throw ("compiling score one strand generates an internal error in Visual C++");
#else
        score_one_strand< DnaAlphabet >(
            pssm,
            BIFA_MAKE_RANGE_FROM_INT_SEQ( sequence ),
            boost::make_function_output_iterator( boost::lambda::bind( &list::append< double >, boost::ref( result ), boost::lambda::_1 ) )
        );
#endif
    } else throw std::logic_error( "Cannot interpret as sequence." );

    return result;
}



/// Functor that takes a log likelihood, position and strandedness and creates a BindingHit
struct make_binding_hit {
    std::string binder_name;
    size_t length;
    biopsy::binding_hit::vec & hits;

    make_binding_hit(
        const std::string & binder_name,
        size_t length,
        biopsy::binding_hit::vec & hits
    )
    : binder_name( binder_name )
    , length( length )
    , hits( hits )
    { }

    void operator()( double log_likelihood, size_t pos, bool positive_strand ) {
        hits.push_back( biopsy::binding_hit( binder_name, biopsy::binding_hit_location( pos, length, positive_strand ), biopsy::get_p_binding( std::exp( log_likelihood ) ) ) );
    }
};


/// Python interface to score both strands
biopsy::binding_hit::vec_ptr
py_score_sequence(
    pssm_t pssm,
    object sequence
) {
    biopsy::binding_hit::vec_ptr result( new biopsy::binding_hit::vec );

    if( can_extract_type< int >( sequence ) ) {
        score_sequence< DnaAlphabet >(
            pssm,
            uniform_sequence_likelihoods(),
            BIFA_MAKE_RANGE_FROM_INT_SEQ( sequence ),
            make_binding_hit( "Binder", boost::size( pssm ), *result )
        );
    } else throw std::logic_error( "Cannot interpret as sequence." );

    return result;
}


} // anonymous namespace




namespace myrrh {
namespace python {
std::string exposed_typechars;
}
}





BOOST_PYTHON_MODULE( _bifa )
{
    using namespace biopsy::bifa;
    using namespace boost::python;

#ifndef NDEBUG
    scope().attr("__debug__") = 1;
    std::cout << "WARNING: Debug version of _bifa module loaded. If you did not intend this then check your configuration!" << std::endl;
#else //_DEBUG
    scope().attr("__debug__") = 0;
#endif //_DEBUG

    scope().attr("__doc__") = "Python interface to C++ library to run BiFa algorithm.";

    import_array();
    myrrh::python::expose_man_fns();
    myrrh::python::expose_converters< double >();

    def(
        "version",
        version,
        "The version of the C++ BiFa python extension" );

    def(
        "score_word",
        py_score_word,
        "Score a PSSM on a word." );

    def(
        "score_one_strand",
        py_score_one_strand,
        "Score a PSSM on the positive strand of a sequence." );

    def(
        "score_sequence",
        py_score_sequence,
        "Score a PSSM on both strands of a sequence." );
}

