/**
@file

Copyright John Reid 2006

*/

#include <boost/python/args.hpp>
#include <boost/python/call.hpp>
#include <boost/python/call_method.hpp>
#include <boost/python/default_call_policies.hpp>
#include <boost/python/suite/indexing/slice_handler.hpp>
#include "biopsy/binding_hits.h"
#include "biopsy/python.h"


using namespace boost;
using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;

#if 0
namespace boost { namespace python { namespace indexing {
template<>
struct value_traits<biopsy::binding_hit> : value_traits< int >
{
	static bool const equality_comparable = false;
	static bool const less_than_comparable = false;
};
} } }
#endif //0

namespace biopsy {


struct hit_location_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs( binding_hit_location const & location )
	{
		return
			boost::python::make_tuple(
				location._position,
				location._length,
				location._positive_strand );
	}
};


struct hit_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs( binding_hit const & hit )
	{
		return
			boost::python::make_tuple(
				hit._binder_name,
				hit._location,
				hit._p_binding );
	}
};

struct hit_vec_pickle_suite : boost::python::pickle_suite
  {
    static
    boost::python::tuple
    getinitargs(binding_hit::vec const& hits)
    {
		return boost::python::tuple();
    }

    static
    boost::python::tuple
    getstate(binding_hit::vec const& hits)
    {
		return boost::python::tuple(hits);
    }

    static
    void
    setstate(binding_hit::vec& hits, boost::python::tuple state)
    {
		hits.clear();
		unsigned len = boost::python::len(state);
		for( unsigned i = 0; len != i; ++i ) hits.push_back( boost::python::extract< binding_hit const& >( state[i] ) );
    }
};

void export_binding_hits()
{
	/**
	A location of a hit
	*/
	class_<
		binding_hit_location
	>(
		"HitLocation",
		"Where a putative hit is located",
		boost::python::init< int, int, bool >( ( boost::python::arg( "position" ), "length", "positive_strand" ) )
	)
		.def_readwrite( "position", &binding_hit_location::_position )
		.def_readwrite( "positive_strand", &binding_hit_location::_positive_strand )
		.def_readwrite( "length", &binding_hit_location::_length )
		.def_pickle( hit_location_pickle_suite() )
		;



	/**
	A hit
	*/
	class_<
		binding_hit
	>(
		"Hit",
		"The details of one putative binding hit",
		boost::python::init< std::string, binding_hit_location, double >()
	)
		.def_readwrite( "binder", &binding_hit::_binder_name )
		.def_readwrite( "location", &binding_hit::_location )
		.def_readwrite( "p_binding", &binding_hit::_p_binding )
		.def_pickle( hit_pickle_suite() )
		;



	/**
	A vector of binding hits.
	*/
	class_<
		binding_hit::vec,
		binding_hit::vec_ptr
	>(
		"HitVec",
		"A sequence of binding hits"
	)
		.def( container_suite< binding_hit::vec >() )
		.def_pickle( hit_vec_pickle_suite() )
		;
	def( "get_binder_names", get_binder_names );



	/**
	A vector of vector of binding hits.
	*/
	class_<
		binding_hits_vec,
		binding_hits_vec_ptr
	>(
		"HitsVec",
		"A sequence of sequences of binding hits"
	)
		.def( container_suite< binding_hits_vec >() )
		;


	/**
	Sort hits.
	*/
	def(
		"sort_hits_by_position",
		sort_hits_by_position,
		"Sorts the hits by position in the sequence" );


	/**
	A mapping from keys to analyses.
	*/
	class_<
		analysis,
		analysis::ptr
	>(
		"Analysis",
		"Maps strings to binding hits"
	)
		.def(
			"get_keys",
			&analysis::get_keys,
			"Gets all the strings with hit sequences associated with them" )
		.def(
			"get_hits_for",
			&analysis::get_hits_for,
			"Gets the hit sequence associated with the key" )
		.def(
			"set_hits_for",
			&analysis::set_hits_for,
			"Sets the hit sequence associated with the key" )
		.def(
			"serialise",
			&analysis::serialise,
			"Serialises to a file" )
		.def(
			"deserialise",
			&analysis::deserialise,
			"Deserialises from a file" )
		.staticmethod( "deserialise" )
		;



}



} //namespace biopsy
