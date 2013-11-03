
#ifndef BIO_PATHWAY_ASSOCIATIONS_H_
#define BIO_PATHWAY_ASSOCIATIONS_H_

#include "bio/common.h"
#include "bio/counter.h"
#include "bio/equivalent_factors.h"
#include "bio/singleton.h"

#include <boost/shared_ptr.hpp>

#include <map>



BIO_NS_START




/** Manages the pathways associated with Matrices, Sites, Factors or Molecules. */
struct PathwayAssociations
	: Singleton< PathwayAssociations >
{
	typedef boost::shared_ptr< PathwayAssociations > ptr_t;
	typedef Counter< TableLink > counter_t;
	typedef std::map< TableLink, counter_t > map_t;
	typedef std::map< EquivalentFactors::partition_ptr_t, counter_t > factor_partition_map_t;

	map_t associations; /**< Holds associations for matrices, sites and individual factors. */
	factor_partition_map_t factor_partition_associations; /**< Holds associations for factor partitions. */

	const counter_t & get_pathways_for(const TableLink & link);
	const counter_t & get_pathways_for(EquivalentFactors::partition_ptr_t factor_partition);
	TableLink get_most_significant_pathway_for(const TableLink & link);

	void init_singleton();
};


BIO_NS_END



#endif //BIO_PATHWAY_ASSOCIATIONS_H_
