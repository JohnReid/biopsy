/**
@file

Copyright John Reid 2006

*/

#include "biopsy/lcs.h"


using namespace boost;
using namespace std;

namespace biopsy {


lcs_ptr lcs_create( const binding_hits_vec & hits_vec )
{
	using namespace boost;

	//create the LCS
	return
		lcs_ptr(
			new lcs(
				make_iterator_range(
					make_indirect_iterator( hits_vec.begin() ),
					make_indirect_iterator( hits_vec.end() ) ) ) );
}

void lcs_calculate( lcs_ptr _lcs )
{
	_lcs->calculate_best();
}


binding_hit::vec_ptr
lcs_get( lcs_ptr _lcs )
{
	//put the result in the vector
	binding_hit::vec_ptr result( new binding_hit::vec );
	BOOST_FOREACH( const lcs::value_holder & v, _lcs->get_best().get_values() )
	{
		result->push_back(
			binding_hit(
				v._c,
				binding_hit_location( v._start, v._end - v._start ),
				v._score ) );
	}

	return result;
}









} //namespace biopsy



