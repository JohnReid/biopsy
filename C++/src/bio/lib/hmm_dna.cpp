/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/environment.h"
#include "bio/hmm_dna.h"
#include "bio/hmm_baum_welch.h"
#include "bio/random.h"
#include "bio/serialisable.h"

#include <boost/iterator/reverse_iterator.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

# ifdef _MSC_VER
# pragma warning( push )
//	warning C4094: untagged 'struct' declared no symbols
# pragma warning( disable : 4094 ) 
# endif
#include <boost/serialization/export.hpp>
# ifdef _MSC_VER
# pragma warning( pop )
# endif
using namespace boost;

#include <vector>
#include <fstream>
using namespace std;

BIO_NS_START

DnaModel::~DnaModel() { }

DnaModel::ptr_t
create_dna_model(unsigned num_states, unsigned order)
{
	if (0 == num_states)
	{
		throw std::invalid_argument("Need positive number of states");
	}

	switch (order)
	{
	case  0: return DnaModel::ptr_t(new DnaHmm< 0>(num_states));
	case  1: return DnaModel::ptr_t(new DnaHmm< 1>(num_states));
#ifndef _DEBUG //hard to debug in MSVC with too many different instantiations
	case  2: return DnaModel::ptr_t(new DnaHmm< 2>(num_states));
	case  3: return DnaModel::ptr_t(new DnaHmm< 3>(num_states));
	case  4: return DnaModel::ptr_t(new DnaHmm< 4>(num_states));
	case  5: return DnaModel::ptr_t(new DnaHmm< 5>(num_states));
	case  6: return DnaModel::ptr_t(new DnaHmm< 6>(num_states));
	case  7: return DnaModel::ptr_t(new DnaHmm< 7>(num_states));
#if 0
	case  8: return DnaModel::ptr_t(new DnaHmm< 8>(num_states));
	case  9: return DnaModel::ptr_t(new DnaHmm< 9>(num_states));
	case 10: return DnaModel::ptr_t(new DnaHmm<10>(num_states));
	case 11: return DnaModel::ptr_t(new DnaHmm<11>(num_states));
	case 12: return DnaModel::ptr_t(new DnaHmm<12>(num_states));
	case 13: return DnaModel::ptr_t(new DnaHmm<13>(num_states));
	case 14: return DnaModel::ptr_t(new DnaHmm<14>(num_states));
	case 15: return DnaModel::ptr_t(new DnaHmm<15>(num_states));
#endif
#endif
	default:
		throw std::logic_error( "Cannot construct DnaHmm with that large an order" );
	}
}





void DnaHmmOrderNumStateMap::init_singleton()
{
	namespace fs = boost::filesystem;

	try_to_deserialise< false >(
		*this,
		fs::path(
			BioEnvironment::singleton().get_species_hmm_file().c_str()
		)
	);
}


BIO_NS_END


#define BIO_DNA_HMM_EXPORT_GUID(O)															\
	typedef BIO_NS::DnaHmm<O> DnaHmm ## O;		\
	BOOST_CLASS_EXPORT(DnaHmm ## O)

//	BOOST_CLASS_EXPORT_GUID(DnaHmm ## O, "DnaHmm" ## BOOST_PP_STRINGIZE(O))

BIO_DNA_HMM_EXPORT_GUID( 0)
BIO_DNA_HMM_EXPORT_GUID( 1)
BIO_DNA_HMM_EXPORT_GUID( 2)
BIO_DNA_HMM_EXPORT_GUID( 3)
BIO_DNA_HMM_EXPORT_GUID( 4)
BIO_DNA_HMM_EXPORT_GUID( 5)
BIO_DNA_HMM_EXPORT_GUID( 6)
BIO_DNA_HMM_EXPORT_GUID( 7)
BIO_DNA_HMM_EXPORT_GUID( 8)
BIO_DNA_HMM_EXPORT_GUID( 9)
BIO_DNA_HMM_EXPORT_GUID(10)
BIO_DNA_HMM_EXPORT_GUID(11)
BIO_DNA_HMM_EXPORT_GUID(12)
BIO_DNA_HMM_EXPORT_GUID(13)
BIO_DNA_HMM_EXPORT_GUID(14)
BIO_DNA_HMM_EXPORT_GUID(15)

