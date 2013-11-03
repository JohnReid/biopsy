
#ifndef BIO_BINDING_MODEL_H_
#define BIO_BINDING_MODEL_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/binding_hit.h"

BIO_NS_START


//forward decl
struct BindingModel;



/**
Allows arbitrary info indexed by type be passed to a BindingModel.
*/
typedef std::map< std::type_info, boost::any > BindingModelContext;












/**
Models something binding to a sequence.
*/
struct BindingModel
{
	typedef boost::shared_ptr< BindingModel > ptr_t;
	typedef BindingHit< BindingModel > hit_t;
	typedef BindingHitSet< BindingModel >::type hit_set_t;
	typedef std::set< BindingModel * > set_t;



	struct parameter_t
	{
		virtual ~parameter_t();
		virtual BindingModel * get_model() const = 0;
	};


	/** Destructor. */
	virtual ~BindingModel();

	/** The name of this model. */
	virtual std::string get_name() const = 0;

	/** The number of bases this model binds to. */
	virtual unsigned get_num_bases() const = 0;

	/** The number of bases this model binds to. */
	virtual const parameter_t * get_parameters() const = 0;

	/** The probability of the binding hypothesis given the sequence. */
	virtual double operator()(
		seq_t::const_iterator begin,
		bool match_complement,
		BindingModelContext * context = 0) const = 0;
};

std::ostream &
operator<<( std::ostream & os, const BindingModel * model );





template< typename is_saving >
struct do_binding_model_serialisation
{
};

//specialise for saving archive
template< >
struct do_binding_model_serialisation< boost::mpl::bool_< true > >
{
	template< typename Archive >
	void operator()( Archive & ar, const BindingModel * & model ) const
	{
		const BindingModel::parameter_t * const parameters = model->get_parameters();
		ar << parameters;
	}

	template< typename Archive >
	void operator()( Archive & ar, BindingModel * & model ) const
	{
		const BindingModel::parameter_t * const parameters = model->get_parameters();
		ar << parameters;
	}
};

//specialise for loading archive
template< >
struct do_binding_model_serialisation< boost::mpl::bool_< false > >
{
	template< typename Archive >
	void operator()( Archive & ar, const BindingModel * & model ) const
	{
		const BindingModel::parameter_t * parameters;
		ar >> parameters;
		model = parameters->get_model();
	}

	template< typename Archive >
	void operator()( Archive & ar, BindingModel * & model ) const
	{
		BindingModel::parameter_t * parameters;
		ar >> parameters;
		model = parameters->get_model();
	}
};


/// Explicit instantiation
template struct do_binding_model_serialisation< boost::mpl::bool_< true > >;

/// Explicit instantiation
template struct do_binding_model_serialisation< boost::mpl::bool_< false > >;







/** Scores a sequence. */
template< typename HitInsIt >
struct SequenceScorer
{
	seq_t::const_iterator seq_begin;
	seq_t::const_iterator seq_end;
	double threshold;
	HitInsIt hit_inserter;
	int start_position;
	BindingModelContext * context;

	SequenceScorer(
		seq_t::const_iterator seq_begin,
		seq_t::const_iterator seq_end,
		double threshold,
		HitInsIt hit_inserter,
		BindingModelContext * context = 0,
		int start_position = 0 )
		: seq_begin( seq_begin )
		, seq_end( seq_end )
		, threshold( threshold )
		, hit_inserter( hit_inserter )
		, start_position( start_position )
		, context( context )
	{
	}

	void operator()(
		BindingModel::ptr_t model )
	{
		return ( *this )( model.get() );
	}

	void operator()(
		BindingModel * model )
	{
		( *this )( model, false );
		( *this )( model, true );
	}

	void operator()(
		BindingModel * model,
		bool complementary)
	{
		//for each position for which there is enough room left to score this model
		const unsigned num_bases = model->get_num_bases();
		int position = start_position;
		for (seq_t::const_iterator s = seq_begin;
			unsigned( seq_end - s ) >= num_bases;
			++s, ++position)
		{
			//what is the prob of binding
			const double p_binding = ( *model )( s, complementary );

			//is it over the threshold
			if ( p_binding > threshold )
			{
				//insert into results
				*hit_inserter++ =
					BindingModel::hit_t(
						model,
						p_binding,
						position,
						num_bases,
						complementary);
			}
		}
	}
};

template< typename HitInsIt >
SequenceScorer< HitInsIt >
make_sequence_scorer(
	seq_t::const_iterator seq_begin,
	seq_t::const_iterator seq_end,
	double threshold,
	HitInsIt hit_inserter,
	BindingModelContext * context = 0,
	int start_position = 0 )
{
	return
		SequenceScorer< HitInsIt >(
			seq_begin,
			seq_end,
			threshold,
			hit_inserter,
			context,
			start_position);
}

BIO_NS_END


template< typename Archive >
Archive &
operator&( Archive & ar, const BIO_NS::BindingModel * & model )
{
	BIO_NS::do_binding_model_serialisation< typename Archive::is_saving >()( ar, model );

	return ar;
}

template< typename Archive >
Archive &
operator&( Archive & ar, BIO_NS::BindingModel * & model )
{
	BIO_NS::do_binding_model_serialisation< typename Archive::is_saving >()( ar, model );

	return ar;
}


#ifdef BOOST_NO_IS_ABSTRACT
BOOST_IS_ABSTRACT( BIO_NS::BindingModel::parameter_t );
#endif //BOOST_NO_IS_ABSTRACT

#endif //BIO_BINDING_MODEL_H_

