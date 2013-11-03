/**
@file

Copyright John Reid 2006, 2008

*/

#include "biopsy/python.h"

namespace biopsy {


void std_exception_translator( std::logic_error const & x )
{
    PyErr_SetString( PyExc_UserWarning, x.what() );
}

void c_string_translator( const char * x )
{
    PyErr_SetString( PyExc_UserWarning, x );
}

void string_translator( const std::string & x )
{
    PyErr_SetString( PyExc_UserWarning, x.c_str() );
}


int major_version() { return 1; }
int minor_version() { return 0; }
std::string version() { return BIOPSY_MAKE_STRING( major_version() << "." << minor_version() ); }

} //namespace biopsy

BOOST_PYTHON_MODULE( _biopsy )
{
	using namespace biopsy;
	using namespace boost::python;

	//register_exception_translator< std::logic_error >( std_exception_translator );
	//register_exception_translator< const char * >( c_string_translator );
	//register_exception_translator< std::string >( string_translator );

#ifndef NDEBUG
	scope().attr("__debug__") = 1;
	std::cout << "WARNING: Debug version of _biopsy module loaded. If you did not intend this then check your configuration!" << std::endl;
#else //_DEBUG
	scope().attr("__debug__") = 0;
#endif //_DEBUG

	scope().attr("__doc__") =
		"Python interface to C++ biopsy library.\r\n"
		"Implements PSSM scanning code (amongst miscellaneous other things).";

	def(
		"version",
		version,
		"The version of the C++ biopsy python extension" );

	export_user();
	export_analyse();
	export_binding_hits();
	export_biopsy();
	// export_db_ref();
	export_lcs();
	export_pssm();
	export_pssm_match();
	export_remo();
	export_test_case();
	export_transfac();
}

