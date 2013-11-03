/**
@file

Copyright John Reid 2006-2009

*/

#include <boost/python.hpp>
#include "biopsy/defs.h"
#include "biopsy/init.h"
#include "biopsy/sequence.h"
#include "biopsy/python.h"

using namespace boost;
using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;



namespace biopsy {


struct sequence_vec_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs( sequence_vec const & o )
	{
		return boost::python::make_tuple();
	}

	static
	boost::python::tuple
	getstate( const sequence_vec & o )
	{
		return boost::python::tuple( o );
	}

	static
	void
	setstate( sequence_vec & o, boost::python::tuple state )
	{
		o.clear();
		unsigned len = boost::python::len(state);
		for( unsigned i = 0; len != i; ++i )
			o.push_back( boost::python::extract< sequence_vec::value_type >( state[i] ) );
	}

};


void export_biopsy()
{
	/**
	A vector of sequences.  The .def (container_suite adds the support of
	*/
	class_< sequence_vec >( "SequenceVec" )
		.def( container_suite< sequence_vec >() )
		.def_pickle( sequence_vec_pickle_suite() )
		;



	/**
	A vector of strings.
	Use register_ptr_to_python so that the shared pointer type is handled corr
	*/
	register_ptr_to_python< string_vec_ptr >();



	/**
	Initialise the module
	*/
	def( "init", biopsy::init );


	/**
	The module version and build
	*/
	def( "get_build", get_build );




	/**
	Reverse complement.
	*/
	def( "reverse_complement", reverse_complement );

	def( "generate_random_sequence", generate_random_sequence );
}



} //namespace biopsy
