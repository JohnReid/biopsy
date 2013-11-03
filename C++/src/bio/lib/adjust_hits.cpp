/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/adjust_hits.h"
#include "bio/binding_model.h"
#include "bio/binding_hit.h"

#include <boost/foreach.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>


BIO_NS_START


struct change_p_binding
{
	change_p_binding( double new_p_binding )
		: new_p_binding( new_p_binding )
	{
	}

	void operator()( BindingHit< BindingModel > & hit ) const
	{
		hit.p_binding = new_p_binding;
	}

private:
	double new_p_binding;
};





void
adjust_hits(
	BindingHitSet< BindingModel >::type & hits,
	seq_t::const_iterator phylo_begin,
	seq_t::const_iterator phylo_end,
	double threshold,
	BindingModelContext * context )
{
	//for each binder in the hits
#ifdef _MSC_VER
	typedef BindingHitSet< BindingModel >::binder binder_t;
	typedef BindingModel::hit_set_t::index< binder_t >::type hit_by_binder_t;
	hit_by_binder_t & by_binder = ::boost::multi_index::get< binder_t >( hits );
#else // _MSC_VER
	typedef BindingHitSet< BindingModel >::type container_t; //should be same as binding_hit_set_t::type
	container_t test; test.begin();
	//BOOST_STATIC_ASSERT( ( boost::is_same<
	//	container_t,
	//	BindingModel::hit_set_t > ) );
	//typedef BindingHitSet< BindingModel > binding_hit_set_t;
	//typedef binding_hit_set_t::binder binder_tag;
	//typedef container_t::index< binder_tag > index_t;
	//typedef index_t::type hit_by_binder_t;
	//hit_by_binder_t & by_binder = ::boost::multi_index::get< binder_tag >( hits );
#endif // _MSC_VER

#if 0
	for( hit_by_binder_t::iterator h = by_binder.begin();
		by_binder.end() != h;
		)
	{
		BindingModel * binder = h->binder;

		//score the model over the phylo sequence
		BindingModel::hit_set_t model_hits;
		make_sequence_scorer(
			phylo_begin,
			phylo_end,
			threshold,
			std::inserter( model_hits, model_hits.begin() ),
			context
		) (
			binder
		);



		//what is the p_binding anywhere in sequence?
		double p_doesnt_bind = 1.0;
		for( BindingModel::hit_set_t::const_iterator mh = model_hits.begin();
			model_hits.end() != mh;
			++mh )
		{
			p_doesnt_bind *= ( 1.0 - mh->p_binding );
		}
		const double p_does_bind = 1.0 - p_doesnt_bind;



		//adjust all the hits for this binder
		while ( by_binder.end() != h && binder == h->binder )
		{
			by_binder.modify( h, change_p_binding( h->p_binding * p_does_bind ) );
			++h;
		}

	}
#endif
}


void
adjust_hits(
	BindingModel::hit_set_t & hits,
	const SeqList & sequences,
	double threshold,
	BindingModelContext * context )
{
	unsigned num_seqs = 0;
	BOOST_FOREACH( const seq_t & s, sequences )
	{
		adjust_hits(
			hits,
			s.begin(),
			s.end(),
			threshold,
			context );
		++num_seqs;
	}

	//adjust all the hits
	for( BindingModel::hit_set_t::const_iterator h = hits.begin();
		hits.end() != h;
		++h )
	{
		hits.modify( 
			h, 
			change_p_binding( 
				exp( log( h->get_p_binding() ) / ( num_seqs + 1 ) ) ) );
	}
}


void
adjust_hits(
	BindingModel::hit_set_t & hits,
	const std::vector< seq_t > & sequences,
	double threshold,
	BindingModelContext * context )
{
	unsigned num_seqs = 0;
	BOOST_FOREACH( const seq_t & s, sequences )
	{
		adjust_hits(
			hits,
			s.begin(),
			s.end(),
			threshold,
			context );
		++num_seqs;
	}

	//adjust all the hits
	for( BindingModel::hit_set_t::const_iterator h = hits.begin();
		hits.end() != h;
		++h )
	{
		hits.modify( 
			h, 
			change_p_binding( 
				exp( log( h->get_p_binding() ) / ( num_seqs + 1 ) ) ) );
	}
}

#if 0
void
remove_under_threshold(
	BindingModel::hit_set_t & hits,
	double threshold )
{
	//for each binder in the hits
	typedef BindingModel::hit_set_t::index< BindingHitSet< BindingModel >::prob >::type hit_by_prob_t;
	hit_by_prob_t & by_prob = hits.get< BindingHitSet< BindingModel >::prob >();

	hit_by_prob_t::const_iterator h = by_prob.begin();
	for( ;
		by_prob.end() != h && h->get_p_binding() < threshold;
		++h )
	{
	}

	by_prob.erase( by_prob.begin(), h );
}
#endif


BIO_NS_END

