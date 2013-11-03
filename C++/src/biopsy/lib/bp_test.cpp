

#include <boost/python.hpp>
#include <boost/python/suite/indexing/container_suite.hpp>

#include <vector>

using namespace boost;
using namespace boost::python;
using namespace std;

typedef vector< string > string_vec;

BOOST_PYTHON_MODULE( _bptest )
{
	/**
	A vector of strings.
	*/
	class_< string_vec >( "StringVec" )
		.def( container_suite< string_vec >() )
		;
}

