/**
Copyright John Reid 2006-2010

@file Code to handle test cases created from TRANSFAC site and fragment tables.

*/

#ifndef BIOPSY_TESTCASE_H_
#define BIOPSY_TESTCASE_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"

namespace biopsy
{

/// A test case. Holds positive and negative sequences.
struct test_case
{
	typedef boost::shared_ptr< test_case > ptr;
	typedef std::vector< ptr > vec;
	typedef boost::shared_ptr< vec > vec_ptr;

	/** Name to identify test case. */
	std::string _name;

	/** Those pssms that should be found in the test case. */
	string_vec_ptr _positives;

	/** Those pssms that should not be found in the test case. */
	string_vec_ptr _negatives;

	/** The sequences in the test case. */
	sequence_vec_ptr _sequences;

	test_case( const std::string & name );
};

test_case::vec_ptr get_transfac_site_test_cases();
test_case::vec_ptr get_transfac_fragment_test_cases();

} //namespace biopsy

#endif //BIOPSY_TESTCASE_H_
