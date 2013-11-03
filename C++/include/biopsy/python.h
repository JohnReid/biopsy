/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_PYTHON_H_
#define BIOPSY_PYTHON_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include <boost/python.hpp>
#include "biopsy/defs.h"

#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/refcount.hpp>
#include <boost/python/suite/indexing/vector.hpp>
#include <boost/python/suite/indexing/list.hpp>


#ifdef _WIN32
#ifdef BIOPSY_BUILD_PYTHON_EXPORT
# define BIOPSY_EXPORT_FN_SPEC __declspec( dllexport )
#else
# define BIOPSY_EXPORT_FN_SPEC __declspec( dllimport )
#endif
#else //_WIN32
# define BIOPSY_EXPORT_FN_SPEC
#endif //_WIN32



namespace biopsy
{

template< typename T >
std::string
as_string( const T & t )
{
	return BIOPSY_MAKE_STRING( t );
}

BIOPSY_EXPORT_FN_SPEC void export_analyse();
BIOPSY_EXPORT_FN_SPEC void export_user();
BIOPSY_EXPORT_FN_SPEC void export_binding_hits();
BIOPSY_EXPORT_FN_SPEC void export_biopsy();
BIOPSY_EXPORT_FN_SPEC void export_db_ref();
BIOPSY_EXPORT_FN_SPEC void export_gapped_pssm();
BIOPSY_EXPORT_FN_SPEC void export_lcs();
BIOPSY_EXPORT_FN_SPEC void export_pssm();
BIOPSY_EXPORT_FN_SPEC void export_pssm_match();
BIOPSY_EXPORT_FN_SPEC void export_remo();
BIOPSY_EXPORT_FN_SPEC void export_test_case();
BIOPSY_EXPORT_FN_SPEC void export_transfac();
BIOPSY_EXPORT_FN_SPEC void export_transfac_2();
BIOPSY_EXPORT_FN_SPEC void export_transfac_3();


void std_exception_translator( std::logic_error const & x );
void c_string_translator( const char * x );
void string_translator( const std::string & x );


inline
boost::python::tuple
boost_to_python_tuple( const boost::tuples::null_type & )
{
	return boost::python::make_tuple();
}

template< typename H, typename T >
inline
boost::python::object
boost_to_python_tuple( const boost::tuples::cons< H, T > & x )
{
	return boost::python::make_tuple( x.get_head() ) + boost_to_python_tuple( x.get_tail() );
}


template< typename T >
struct tupleconverter
{
	static PyObject * convert( T const & x )
	{
		return boost::python::incref( boost_to_python_tuple( x ).ptr() );
	}
};



} //namespace biopsy

#endif //BIOPSY_PYTHON_H_


