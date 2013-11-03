/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/bifa_analysis.h"
#include "bio/biobase_binding_model.h"



BIO_NS_START


void bifa_hits_2_match_results(
	const bifa_hits_t & hits,
	match_result_vec_t & results,
	double threshold )
{
	BOOST_FOREACH( const BindingHit< BindingModel > & hit, ::boost::multi_index::get< BindingHitSet< BindingModel >::position >( hits ) )
	{
		if( hit.p_binding > threshold )
		{
			const BiobaseBindingModel * binding_model = dynamic_cast< const BiobaseBindingModel * >( hit.binder );
			if( 0 != binding_model )
			{
				match_result_vec_t::value_type result(
					binding_model->parameters.link,
					Hit(
						float_t( hit.p_binding ),
						hit.position,
						hit.complementary ) );

				results.push_back( result );
			}
		}
	}
}


BIO_NS_END

