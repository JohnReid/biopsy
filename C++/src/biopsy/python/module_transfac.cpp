/**
@file

Copyright John Reid 2006, 2008

*/

#include "biopsy/python.h"

namespace biopsy {



int major_version() { return 1; }
int minor_version() { return 0; }
std::string version() { return BIOPSY_MAKE_STRING( major_version() << "." << minor_version() ); }

} //namespace biopsy

BOOST_PYTHON_MODULE( _transfac )
{
	using namespace biopsy;
	using namespace boost::python;

#ifndef NDEBUG
	scope().attr("__debug__") = 1;
	std::cout << "WARNING: Debug version of _transfac module loaded. If you did not intend this then check your configuration!" << std::endl;
#else //_DEBUG
	scope().attr("__debug__") = 0;
#endif //_DEBUG

	scope().attr("__doc__") =
		"Python interface to C++ library to access TRANSFAC data.";

	def(
		"version",
		version,
		"The version of the C++ TRANSFAC python extension" );

	export_transfac_2();
}

