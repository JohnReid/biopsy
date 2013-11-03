
#ifndef BIO_MATCH_BINDING_MODEL_H_
#define BIO_MATCH_BINDING_MODEL_H_

#include "bio/defs.h"
#include "bio/environment.h"
#include "bio/common.h"
#include "bio/biobase_pssm.h"
#include "bio/binding_model.h"

#include <boost/iterator/transform_iterator.hpp>
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


/** The match algorithm. */
struct MatchBindingModel
	: BindingModel
{
	struct parameter_t
		: BindingModel::parameter_t
		, boost::less_than_comparable< MatchBindingModel::parameter_t >
	{
		TableLink link;

		parameter_t(
			TableLink link = TableLink() );

		virtual ~parameter_t();
		
		virtual BindingModel * get_model() const;
		
		bool operator<( const MatchBindingModel::parameter_t & rhs ) const;

	private:
		friend class boost::serialization::access;
		template< typename Archive >
		void serialize( Archive & ar, const unsigned int version )
		{
			boost::serialization::void_cast_register< MatchBindingModel::parameter_t, BindingModel::parameter_t >( 0, 0 );

			ar & link;
		}
	};

	parameter_t parameters;
	Pssm pssm;
	double threshold;

	MatchBindingModel( const MatchBindingModel::parameter_t & parameters = parameter_t() );

	/** Destructor. */
	virtual ~MatchBindingModel();

	/** The name of this model. */
	virtual std::string get_name() const;

	/** The number of bases this model binds to. */
	virtual unsigned get_num_bases() const;

	/** The parameters for this model. */
	virtual const parameter_t * get_parameters() const;

	/** The probability of the binding hypothesis given the sequence. */
	virtual double operator()(
		seq_t::const_iterator begin,
		bool match_complement,
		BindingModelContext * context) const;
};


struct Link2MatchBindingModel
{
	BindingModel * operator()( const TableLink & link ) const
	{
		return MatchBindingModel::parameter_t( link ).get_model();
	}

	template< typename BiobaseMapValue >
	BindingModel * operator()( const BiobaseMapValue & value ) const
	{
		return ( *this )( value.first );
	}
};




BIO_NS_END

#endif //BIO_MATCH_BINDING_MODEL_H_
