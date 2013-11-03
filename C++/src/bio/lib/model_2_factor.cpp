/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/model_2_factor.h"
#include "bio/transcription_factor.h"
#include "bio/biobase_binding_model.h"
#include "bio/biobase_likelihoods.h"
#include "bio/biobase_tf.h"
#include "bio/biobase_db.h"

#include <boost/function_output_iterator.hpp>

BIO_NS_START




template< typename OutputIt >
void
get_factors_from_model( const BindingModel * model, OutputIt it )
{
	const BiobaseBindingModel * biobase_binding_model = dynamic_cast< const BiobaseBindingModel * >( model );
	if ( 0 != biobase_binding_model )
	{
		BiobaseTablePssmEntry * pssm = BiobaseDb::singleton().get_pssm_entry( biobase_binding_model->parameters.link );
		EquivalentFactors::partition_set_t factors = EquivalentFactors::singleton().get_factors_for( pssm );

		std::transform(
			factors.begin(),
			factors.end(),
			it,
			boost::bind< TF & >( boost::ref( BiobaseTFCache::singleton() ), _1 )
		);
	}
}


template< typename OutputIt >
struct ModelHit2FactorHits
	: std::unary_function< BindingModel::hit_t, void >
{
	OutputIt output_it;

	ModelHit2FactorHits(
		OutputIt output_it )
		: output_it( output_it )
	{
	}

	result_type operator()( argument_type model_hit )
	{
		get_factors_from_model( 
			model_hit.binder,
			boost::make_function_output_iterator( FactorHitCreator( model_hit, output_it ) ) );
	}

	struct FactorHitCreator
	{
		const BindingModel::hit_t & model_hit;
		OutputIt output_it;

		FactorHitCreator(
			const BindingModel::hit_t & model_hit,
			OutputIt output_it )
			: model_hit( model_hit )
			, output_it( output_it )
		{
		}

		void operator()( const TF & factor )
		{
			*output_it =
				TF::hit_t(
					const_cast< TF * >( boost::addressof( factor ) ),
					model_hit.p_binding,
					model_hit.position,
					model_hit.length,
					model_hit.complementary );

			output_it++;
		}
	};
};

template< typename OutputIt >
ModelHit2FactorHits< OutputIt >
make_model_2_factors( OutputIt output_it )
{
	return ModelHit2FactorHits< OutputIt >( output_it );
}




template< typename Binder, typename OutputIt >
void
collapse_hits(
	typename BindingHitSet< Binder >::type::template index< typename BindingHitSet< Binder >::binder >::type::const_iterator begin,
	typename BindingHitSet< Binder >::type::template index< typename BindingHitSet< Binder >::binder >::type::const_iterator end,
	OutputIt output_it )
{
	typedef BindingHit< Binder > hit_t;
	typedef typename BindingHitSet< Binder >::type hit_set_t;
	typedef typename hit_set_t::template index< typename BindingHitSet< Binder >::binder >::type hit_set_by_binder_t;
	typedef typename hit_set_by_binder_t::const_iterator const_iterator;

	//iterate through by position
	while( end != begin )
	{
		hit_t hit = *begin;
		++begin;
		int hit_end = hit.get_end();

		//while we are still on the same binder and overlapping keep extending the end
		//and updating the binding prob and the complementarity
		while ( end != begin
			&& hit.get_end() > begin->get_position() 
			&& hit.get_binder() == begin->get_binder() )
		{
			hit_end = std::max( hit_end, begin->get_end() );
			if ( hit.p_binding < begin->get_p_binding() )
			{
				hit.p_binding = begin->get_p_binding();
				hit.complementary = begin->is_complementary();
			}
			++begin;
		}
		hit.length = hit_end - hit.get_position();

		*output_it = hit;
		output_it++;
	}
}


void model_hits_2_factor_hits(
	const BindingModel::hit_set_t & model_hits,
	TF::hit_set_t & factor_hits )
{
	TF::hit_set_t uncollapsed_hits;
	std::for_each(
		model_hits.begin(),
		model_hits.end(),
		make_model_2_factors( std::inserter( uncollapsed_hits, uncollapsed_hits.begin() ) ) );

	factor_hits.clear();
	collapse_hits< TF >(
		uncollapsed_hits.get< BindingHitSet< TF >::binder >().begin(),
		uncollapsed_hits.get< BindingHitSet< TF >::binder >().end(),
		std::inserter( factor_hits, factor_hits.begin() ) );
}



BIO_NS_END




