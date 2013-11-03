/**
@file

Copyright John Reid 2006

*/

#include "biopsy/python.h"
#include "biopsy/lcs.h"

using namespace boost;
using namespace boost::python;
using namespace std;


namespace biopsy {

void export_lcs()
{
	/**
	Longest common subsequence.
	*/
	class_< lcs >( "LCS", no_init )
		.def( "get_storage_size", &lcs::get_storage_size )
		;
	register_ptr_to_python< lcs_ptr >();

	def( "lcs_create", lcs_create );
	def( "lcs_calculate", lcs_calculate );
	def( "lcs_get", lcs_get );



}



} //namespace biopsy
