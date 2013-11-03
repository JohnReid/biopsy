
#ifndef BIO_BIOBASE_SCORE_H_
#define BIO_BIOBASE_SCORE_H_

#include "bio/defs.h"
#include "bio/environment.h"
#include "bio/common.h"
#include "bio/biobase_filter.h"
#include "bio/biobase_pssm.h"
#include "bio/biobase_binding_model.h"
#include "bio/biobase_db.h"

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/function_output_iterator.hpp>

BIO_NS_START




template<
	typename Filter,
	typename Transformer,
	typename OutputIt
>
void
transform_biobase_sites_and_matrices(
	Filter filter,
	Transformer transformer,
	OutputIt output_it )
{
	using namespace boost;

	std::transform(
		make_filter_iterator(
			filter,
			BiobaseDb::singleton().get_matrices().begin(),
			BiobaseDb::singleton().get_matrices().end()
		),
		make_filter_iterator(
			filter,
			BiobaseDb::singleton().get_matrices().end(),
			BiobaseDb::singleton().get_matrices().end()
		),
		output_it,
		transformer
	);

	std::transform(
		make_filter_iterator(
			filter,
			BiobaseDb::singleton().get_sites().begin(),
			BiobaseDb::singleton().get_sites().end()
		),
		make_filter_iterator(
			filter,
			BiobaseDb::singleton().get_sites().end(),
			BiobaseDb::singleton().get_sites().end()
		),
		output_it,
		transformer
	);
}



/**
Scores all the biobase pssms that match the filter.
*/
template<
	typename SeqScorer,
	typename Link2BindingModel,
	typename LinkFilter
>
void
score_all_biobase_pssms(
	const SeqScorer & sequence_scorer,
	LinkFilter filter = BiobasePssmFilter(),
	Link2BindingModel link_2_binding_model = Link2BiobaseBindingModel() )
{
	transform_biobase_sites_and_matrices(
		filter,
		link_2_binding_model,
		boost::make_function_output_iterator( sequence_scorer ) );
}

BIO_NS_END

#endif //BIO_BIOBASE_SCORE_H_
