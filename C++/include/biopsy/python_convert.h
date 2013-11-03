/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_PYTHON_CONVERT_H_
#define BIOPSY_PYTHON_CONVERT_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/python.h"
#include "biopsy/numpy_converter.h"
#include "biopsy/containers.h"

#include <boost/multi_array.hpp>

#include <vector>





namespace biopsy {

typedef std::vector< double > double_vector; /**< A vector of double. */
typedef std::vector< double_vector > double_vector_vec; /**< A vector of vector of floats. */
typedef boost::multi_array< double, 2 > double_array;

extern numpy_converter converter;


//forward decl
/** Convert given C++ Class boost::python::object, o, to python boost::python::object. */
template< typename Class >
boost::python::object
convert_to_python( const Class & o );

//forward decl
template< typename Class >
void
convert_from_python( boost::python::object py_obj, Class & o );

namespace impl {

template< typename Class >
struct _convert_to_python
{
	boost::python::object operator()( const Class & o ) const
	{
		BOOST_STATIC_ASSERT( ! "_convert_to_python not specialised for this class" );
		return boost::python::object();
	}
};

template< >
struct _convert_to_python< double_array >
{
	boost::python::object operator()( const double_array & o ) const
	{
		boost::python::object result = converter.to_numpy( o );
		result.attr( "setflags" )( false ); //make not writeable
		return result;
	}
};

template< typename T >
struct _convert_to_python< std::vector< T > >
{
	boost::python::object operator()( const std::vector< T > & o ) const
	{
		boost::python::object result = converter.to_numpy( o );
		result.attr( "setflags" )( false ); //make not writeable
		return result;
	}
};

template< unsigned N >
struct _convert_to_python< boost::array< double, N > >
{
	boost::python::object operator()( const boost::array< double, N > & o ) const
	{
		boost::python::object result = converter.to_numpy( o );
		result.attr( "setflags" )( false ); //make not writeable
		return result;
	}
};

template< typename T >
struct _convert_to_python< std::vector< std::vector< T > > >
{
	boost::python::object operator()( const std::vector< std::vector< T > > & o ) const
	{
		boost::python::list l;
		for( unsigned i = 0; o.size() != i; ++i )
		{
			l.append( convert_to_python( o[i] ) );
		}
		return l;
	}
};

template< >
struct _convert_to_python< double_vector_vec_vec >
{
	boost::python::object operator()( const double_vector_vec_vec & o ) const
	{
		boost::python::list l;
		for( unsigned i = 0; o.size() != i; ++i )
		{
			l.append( convert_to_python( o[i] ) );
		}
		return l;
	}
};

template< typename Class >
struct _convert_from_python
{
	void operator()( boost::python::object py_obj, Class & o ) const
	{
		BOOST_STATIC_ASSERT( ! "_convert_from_python not specialised for this class" );
	}
};

template< >
struct _convert_from_python< double_array >
{
	void operator()( boost::python::object py_obj, double_array & o ) const
	{
		converter.from_numpy( py_obj, o );
	}
};

template< typename T >
struct _convert_from_python< std::vector< T > >
{
	void operator()( boost::python::object py_obj, std::vector< T > & o ) const
	{
		converter.from_numpy( py_obj, o );
	}
};

template< typename T, unsigned N >
struct _convert_from_python< boost::array< T, N > >
{
	void operator()( boost::python::object py_obj, boost::array< T, N > & o ) const
	{
		converter.from_numpy( py_obj, o );
	}
};

template< typename T >
struct _convert_from_python< std::vector< std::vector< T > > >
{
	void operator()( boost::python::object py_obj, std::vector< std::vector< T > > & o ) const
	{
		o.resize( boost::python::len( py_obj ) );
		for( unsigned i = 0; o.size() != i; ++i )
		{
			convert_from_python( py_obj[ i ], o[i] );
		}
	}
};

template< typename T >
struct _convert_from_python< std::vector< std::vector< std::vector< T > > > >
{
	void operator()( boost::python::object py_obj, std::vector< std::vector< std::vector< T > > > & o ) const
	{
		o.resize( boost::python::len( py_obj ) );
		for( unsigned i = 0; o.size() != i; ++i )
		{
			convert_from_python( py_obj[ i ], o[i] );
		}
	}
};


} //namespace impl



template< typename Class >
boost::python::object
convert_to_python( const Class & o )
{
	return impl::_convert_to_python< Class >()( o );
}

template< 
	typename Class,
	typename Container,
	Container Class:: * Member
>
boost::python::object
access_and_convert( const Class & obj )
{
	return convert_to_python( obj.*Member );
}

template< typename Class >
void
convert_from_python( boost::python::object py_obj, Class & o )
{
	impl::_convert_from_python< Class >()( py_obj, o );
}

template< 
	typename Class,
	typename Container,
	Container Class:: * Member
>
void
convert_and_set( Class & obj, boost::python::object py_obj )
{
	convert_from_python( py_obj, obj.*Member );
}



} //namespace biopsy

#endif //BIOPSY_PYTHON_CONVERT_H_
