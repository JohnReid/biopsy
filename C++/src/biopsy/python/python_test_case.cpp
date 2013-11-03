/**
@file

Copyright John Reid 2006

*/

#include <boost/python.hpp>
#include "biopsy/test_case.h"
#include "biopsy/python.h"


using namespace boost;
using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;

namespace boost { namespace python { namespace indexing {
template<>
struct value_traits<biopsy::test_case::ptr> : indirect_value_traits<biopsy::test_case::ptr>
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
} } }



namespace biopsy {



void export_test_case()
{
	class_<
		test_case,
		test_case::ptr
	>(
		"TestCase",
		"A test case for pssms",
		no_init
	)
		.def_readwrite( "name", &test_case::_name )
		.def_readwrite( "positives", &test_case::_positives )
		.def_readwrite( "negatives", &test_case::_negatives )
		.def_readwrite( "sequences", &test_case::_sequences )
		;

	class_<
		test_case::vec,
		test_case::vec_ptr
	>(
		"TestCaseVec",
		"A sequence of pssm test cases"
	)
		.def( container_suite< test_case::vec >() )
		;

	def( "get_transfac_site_test_cases", get_transfac_site_test_cases );
	def( "get_transfac_fragment_test_cases", get_transfac_fragment_test_cases );
}




} //namespace biopsy
