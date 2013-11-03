/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/biobase_tf.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"

# ifdef _MSC_VER
# pragma warning( push )
//	warning C4094: untagged 'struct' declared no symbols
# pragma warning( disable : 4094 ) 
# endif
#include <boost/serialization/export.hpp>
# ifdef _MSC_VER
# pragma warning( pop )
# endif

BIO_NS_START



/**
A transcription factor from biobase.
*/
struct BiobaseTF
	: TF
{
	EquivalentFactors::partition_ptr_t factor_partition;

	BiobaseTF( const EquivalentFactors::partition_ptr_t & factor_partition = EquivalentFactors::partition_ptr_t() );

	virtual std::string get_name() const;

private:
    friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int version) const
	{
		const unsigned indicative_acc_id = EquivalentFactors::singleton().get_indicative_acc_id( factor_partition );
		ar & indicative_acc_id;
	}
	template<class Archive>
	void load(Archive & ar, const unsigned int version)
	{
		unsigned indicative_acc_id;
		ar & indicative_acc_id;
		factor_partition = EquivalentFactors::singleton().get_partition( indicative_acc_id );
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};



BiobaseTF::BiobaseTF( const EquivalentFactors::partition_ptr_t & factor_partition )
: factor_partition( factor_partition )
{
}


std::string
BiobaseTF::get_name() const
{
	if( 0 == factor_partition )
	{
		throw std::logic_error( "Null pointer in BiobaseTF::get_name()" );
	}
	return EquivalentFactors::get_name_for( factor_partition );
}


BiobaseFactor2TF::result_type
BiobaseFactor2TF::operator()( argument_type factor ) const
{
	return BiobaseFactor2TF::result_type( new BiobaseTF( factor ) );
}




BIO_NS_END

BOOST_CLASS_EXPORT( BIO_NS::BiobaseTF )

