/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/pathway_associations.h"
#include "bio/biobase_db.h"
#include "bio/biobase_filter.h"
#include "bio/biobase_data_traits.h"
#include "bio/log.h"
USING_BIO_NS;

using namespace boost;


#include <iostream>
using namespace std;





template <class Functor>
bool
apply_to_interesting_ancestor_pathways(
	Pathway * pathway,
	Functor & functor)
{
	bool found_interesting_pathway = false;

	if (PATHWAY_PW == pathway->pathway_type && IsInterestingPathway()(pathway->accession_number))
	{
		//found interesting pathway
		functor(pathway->accession_number);
		found_interesting_pathway = true;
	}
	else
	{
		//look for interesting ancestor
		for (typename TableLinkVec::const_iterator sf = pathway->super_families.begin();
			pathway->super_families.end() != sf;
			++sf)
		{
			found_interesting_pathway |=
				apply_to_interesting_ancestor_pathways(
					BiobaseDb::singleton().get_entry<PATHWAY_DATA>(*sf),
					functor);
		}
	}

	return found_interesting_pathway;
}


/** Update the counts molecules, factor partitions and factors. */
struct AssociationUpdater
{
	PathwayAssociations & associations;
	PathwayAssociations::counter_t & counter;

	AssociationUpdater(PathwayAssociations & associations, PathwayAssociations::counter_t & counter)
		: associations(associations)
		, counter(counter)
	{
	}

	void operator()(TableLink pathway_link)
	{
		BOOST_ASSERT(PATHWAY_DATA == pathway_link.table_id);

		Pathway * pathway = BiobaseDb::singleton().get_entry<PATHWAY_DATA>(pathway_link);

		if (0 != pathway)
		{
			apply_to_interesting_ancestor_pathways(pathway, counter);
		}
	}

	void operator()(DatabaseRef db_ref)
	{
		//is it a link to TransPath?
		if (TRANSPATH_DB == db_ref.db)
		{
			TableLink molecule_link = transfac_table_link_from_db_ref(db_ref);
			if (MOLECULE_DATA == molecule_link.table_id)
			{
				counter += associations.get_pathways_for(molecule_link);
			}
		}
	}

	void operator()(EquivalentFactors::partition_ptr_t partition)
	{
		counter += associations.get_pathways_for(partition);
	}

	void operator()(unsigned factor_accession_number)
	{
		counter += associations.get_pathways_for(TableLink( FACTOR_DATA, factor_accession_number ));
	}
};



const PathwayAssociations::counter_t &
PathwayAssociations::get_pathways_for(EquivalentFactors::partition_ptr_t factor_partition)
{
	factor_partition_map_t::iterator i = factor_partition_associations.find(factor_partition);
	if (factor_partition_associations.end() == i)
	{
		i = factor_partition_associations.insert(factor_partition_map_t::value_type(factor_partition, counter_t())).first;

		std::for_each(factor_partition->begin(), factor_partition->end(), AssociationUpdater(*this, i->second));
	}
	return i->second;
}


const PathwayAssociations::counter_t &
PathwayAssociations::get_pathways_for(const TableLink & link)
{
	/**
	For each Table link a cache is maintained of the pathways to which it is associated.
	There is potentially a single entry in the assocations map for each matrix, site, factor or molecule
	entry.  If it does not exist then it is created.
		For transcription factors, we find the list of pathways they are associated with
		For Pssms and sites we find the associated factors, and then find the associated pathways
		For molecules...
	*/
		
	map_t::iterator i = associations.find(link);
	if (associations.end() == i)
	{
		i = associations.insert(map_t::value_type(link, counter_t())).first;

		switch(link.table_id)
		{
		case MATRIX_DATA:
		case SITE_DATA:
			{
				BiobaseTablePssmEntry * pssm = BiobaseDb::singleton().get_pssm_entry(link);
				//for each factor
				EquivalentFactors::partition_set_t partitions = EquivalentFactors::singleton().get_factors_for(pssm);
				std::for_each(partitions.begin(), partitions.end(), AssociationUpdater(*this, i->second));
			}
			break;

		case FACTOR_DATA:
			{
				Factor * factor = BiobaseDb::singleton().get_entry<FACTOR_DATA>(link);
				if (0 != factor)
				{
					std::for_each(factor->database_refs.begin(), factor->database_refs.end(), AssociationUpdater(*this, i->second));
				}
			}
			break;

//		case MOLECULE_DATA:
//			{
//				Molecule * molecule = BiobaseDb::singleton().get_entry<MOLECULE_DATA>(link);
//				if (0 != molecule)
//				{
//					std::for_each(molecule->pathways.begin(), molecule->pathways.end(), AssociationUpdater(*this, i->second));
//
//					//if we didn't find any pathways for this molecule look in the super families
//					if (0 == i->second.get_total())
//					{
//						for (TableLinkVec::const_iterator sf = molecule->super_families.begin();
//							molecule->super_families.end() != sf;
//							++sf)
//						{
//							i->second += get_pathways_for(*sf);
//						}
//					}
//				} //if 0 != molecule
//			}
//			break;

		default:
			throw std::logic_error( "Cannot get pathways for this type of data" );
		}
	}
	return i->second;
}

struct PathwayAssociationBuilder
{
	PathwayAssociations & associations;
	unsigned num_with_pathways;
	unsigned total;
	std::string type;

	PathwayAssociationBuilder(PathwayAssociations & associations, std::string type)
		: associations(associations)
		, num_with_pathways(0)
		, total(0)
		, type(type)
	{
	}

	template < typename T >
	void operator()( const T & value )
	{
		++total;
		if (! associations.get_pathways_for(value.second->accession_number).empty())
		{
			++num_with_pathways;
		}
	}

	void print() const
	{
		log_stream() << "# " << type << " with pathways = " << num_with_pathways << " out of " << total << "\n";
	}
};

void PathwayAssociations::init_singleton()
{
	std::for_each(
		get_matrices_begin(),
		get_matrices_end(),
		PathwayAssociationBuilder(*this, "matrices")).print();

	std::for_each(
		get_sites_begin(),
		get_sites_end(),
		PathwayAssociationBuilder(*this, "sites")).print();
}

/** Find the link to the most significant pathway for the given link. Returns UNKNOWN_DATA if none found. */
TableLink
PathwayAssociations::get_most_significant_pathway_for(const TableLink & link)
{
	TableLink result;

	counter_t pathway_hits = get_pathways_for(link);
	counter_t::const_iterator hit_it = pathway_hits.find_max();
	if (pathway_hits.end() != hit_it)
	{
		result = hit_it->first;
	}

	return result;
}





