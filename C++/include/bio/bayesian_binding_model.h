
#ifndef BIO_BAYESIAN_BINDING_MODEL_H_
#define BIO_BAYESIAN_BINDING_MODEL_H_

#include "bio/defs.h"
#include "bio/binding_model.h"

#include <boost/numeric/interval.hpp>


BIO_NS_START




/**
A binding model that relies on the evaluation of some statistic (or data) e.g. the biobase score of a sequence.

Seq2Data::Seq2Data must be convertible to BackgroundLikelihood::first_argument_type and BindingLikelihood::first_argument_type

BackgroundLikelihood::result_type must be convertible to double

BindingLikelihood::result_type must be convertible to double

*/
template<
	typename Seq2Data,				//converts sequences to data
	typename BackgroundLikelihood,	//evaluates the likelihood of the data under the background hypothesis
	typename BindingLikelihood		//evaluates the likelihood of the data under the binding hypothesis
>
struct BayesianBindingModel
	: BindingModel
{
	typedef typename Seq2Data::result_type data_type;

	std::string name;
	Seq2Data seq_2_data;
	BackgroundLikelihood background_likelihood;
	BindingLikelihood binding_likelihood;
	double p_Hb_prior;

	BayesianBindingModel(
		const std::string & name = "uninitialised",
		const Seq2Data & seq_2_data = Seq2Data(), 
		const BackgroundLikelihood & background_likelihood = BackgroundLikelihood(),
		const BindingLikelihood & binding_likelihood = BindingLikelihood(),
		double p_Hb_prior = 0.5)
		: name( name )
		, seq_2_data( seq_2_data )
		, background_likelihood( background_likelihood )
		, binding_likelihood( binding_likelihood )
		, p_Hb_prior( p_Hb_prior )
	{
		using namespace boost::numeric;
		static interval< double > interval_0_1( 0.0, 1.0 ); 

		if ( ! in( p_Hb_prior, interval_0_1 ) ) //should be in [0,1]
		{
			throw std::logic_error( BIO_MAKE_STRING( get_name() << ": p_Hb_prior is out of range" ) );
		}
	}

	/** The number of bases this model binds to. */
	virtual unsigned get_num_bases() const
	{
		return seq_2_data.get_num_bases();
	}

	virtual std::string get_name() const
	{
		return name;
	}

	virtual const parameter_t * get_parameters() const
	{
		throw std::logic_error( "No parameters for this type of model" );
	}

	/** The probability of the binding hypothesis given the sequence. */
	virtual double operator()(
		seq_t::const_iterator begin,
		bool match_complement,
		BindingModelContext * context = 0 ) const
	{
		using namespace boost::numeric;
		static interval< double > interval_0_1( 0.0, 1.0 ); 

		//the data generated from the sequence
		data_type data = seq_2_data( begin, match_complement );

		const double p_D_given_Hb = binding_likelihood( data );
		if ( ! in( p_D_given_Hb, interval_0_1 ) ) //should be in [0,1]
		{
			throw std::logic_error( BIO_MAKE_STRING( get_name() << ": p_D_given_Hb is out of range" ) );
		}

		//this could be 0
		const double p_D_given_Hn = background_likelihood( data );
		if ( ! in( p_D_given_Hn, interval_0_1 ) )
		{
			throw std::logic_error( BIO_MAKE_STRING( get_name() << ": p_D_given_Hn is out of range" ) );
		}
		
		//so this could be infinite - we set to 0 if p_D_given_Hb = 0
		const float_t bayes_factor = float_t( 0 == p_D_given_Hb ? 0.0 : (p_D_given_Hb * p_Hb_prior) / (p_D_given_Hn * ( 1.0 - p_Hb_prior ) ) );

		//in which case this should be 1
		const double p_Hb_given_D =
			BIO_FINITE( bayes_factor )
				? bayes_factor / (1.0 + bayes_factor)
				: 1.0;

		//our prob should be between 0 and 1
		if ( ! in( p_Hb_given_D, interval_0_1 ) )
		{
			throw std::logic_error( BIO_MAKE_STRING( get_name() << ": p_Hb_given_D is out of range" ) );
		}

		return p_Hb_given_D;
	}
};





BIO_NS_END

#endif //BIO_BAYESIAN_BINDING_MODEL_H_
