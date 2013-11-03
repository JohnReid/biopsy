
#ifndef BIO_BIOBASE_BINDING_MODEL_H_
#define BIO_BIOBASE_BINDING_MODEL_H_

#include "bio/defs.h"
#include "bio/environment.h"
#include "bio/common.h"
#include "bio/biobase_filter.h"
#include "bio/biobase_pssm.h"
#include "bio/binding_model.h"
#include "bio/bayesian_binding_model.h"
#include "bio/biobase_likelihoods.h"

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





struct BiobaseBindingModel
	: BayesianBindingModel< PssmScorer, QuantisedScores, QuantisedScores >
{
	struct parameter_t
		: BindingModel::parameter_t
		, boost::less_than_comparable< BiobaseBindingModel::parameter_t >
	{
		TableLink link;
		double p_Hb_prior;
		bool or_better;

		parameter_t(
			const TableLink & link = TableLink(),
			double p_Hb_prior = BioEnvironment::singleton().get_tf_binding_prior(),
			bool or_better = false );

		virtual BindingModel * get_model() const;

		bool operator<( const BiobaseBindingModel::parameter_t & rhs ) const;

	private:
		friend class boost::serialization::access;
		template< typename Archive >
		void serialize( Archive & ar, const unsigned int version )
		{
			boost::serialization::void_cast_register< BiobaseBindingModel::parameter_t, BindingModel::parameter_t >( 0, 0 );

			ar & link;
			ar & p_Hb_prior;
			ar & or_better;
		}
	};


	parameter_t parameters;

	BiobaseBindingModel( const BiobaseBindingModel::parameter_t & parameters = parameter_t() );

	/** The parameters for this model. */
	virtual const BindingModel::parameter_t * get_parameters() const;

private:
	BindingModel & get_underlying_model() const;
};


struct Link2BiobaseBindingModel
{
	double p_Hb_prior;
	bool or_better;

	Link2BiobaseBindingModel( 
		double p_Hb_prior = BioEnvironment::singleton().get_tf_binding_prior(),
		bool or_better = false )
		: p_Hb_prior( p_Hb_prior )
		, or_better( or_better )
	{
	}

	BindingModel * operator()( const TableLink & link ) const
	{
		return BiobaseBindingModel::parameter_t( link, p_Hb_prior, or_better ).get_model();
	}

	template< typename BiobaseMapValue >
	BindingModel * operator()( const BiobaseMapValue & value ) const
	{
		return ( *this )( value.first );
	}
};


BIO_NS_END

#endif //BIO_BIOBASE_BINDING_MODEL_H_
