/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/equivalent_factors.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/environment.h"
#include "bio/serialisable.h"

#include <boost/filesystem/fstream.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/progress.hpp>
namespace fs = boost::filesystem;

#include <iostream>


BIO_NS_START


void
get_all_pssms(pssm_set & pssms)
{
	//
	// Construct a set of all the pssms we know about...
	//

	//put each matrix in the pssm set
	for (Matrix::map_t::const_iterator m = BiobaseDb::singleton().get_matrices().begin();
		BiobaseDb::singleton().get_matrices().end() != m;
		++m)
	{
		pssms.insert(m->second.get());
	}

	//put each consensus site in the pssm set
	for (Site::map_t::const_iterator s = BiobaseDb::singleton().get_sites().begin();
		BiobaseDb::singleton().get_sites().end() != s;
		++s)
	{
		if ("CONS" == s->second->id.factor)
		{
			pssms.insert(s->second.get());
		}
	}

	std::cout << "Found " << pssms.size() << " pssms\n";
}




Factor *
FactorEquivalence::get_factor(unsigned factor_acc_number)
{
	return BiobaseDb::singleton().get_entry<FACTOR_DATA>(TableLink(FACTOR_DATA, factor_acc_number));
}




EquivalentFactors::ptr_t EquivalentFactors::construct_from_biobase()
{
	EquivalentFactors::ptr_t result(new EquivalentFactors);

	result->init_from_biobase();

	return result;
}




void EquivalentFactors::init_from_biobase()
{
	std::cout << "Building list of factor synonyms from Biobase\n";
	boost::progress_timer timer;

	//remove any existing partitions
	partitions.clear();
	
	//
	// Construct a set of all the pssms we know about...
	//
	pssm_set pssms;
	get_all_pssms(pssms);


	//
	// Look at each pssm
	//
	for (pssm_set::const_iterator p = pssms.begin();
		pssms.end() != p;
		++p)
	{
		const FactorLinkList & factors = (*p)->get_factors();

		//for each factor
		for (FactorLinkList::const_iterator f = factors.begin();
			factors.end() != f;
			++f)
		{
			//if it exists in biobase
			if (0 != BiobaseDb::singleton().get_entry<FACTOR_DATA>(f->get()->link))
			{
				add_object(f->get()->link.entry_idx);
			}
		}
	}
}

unsigned
EquivalentFactors::get_indicative_acc_id( partition_ptr_t partition )
{
	if( 0 == partition )
	{
		throw std::logic_error( "Null pointer in EquivalentFactors::get_indicative_acc_id()" );
	}

	return *( partition->begin() );
}


void
EquivalentFactors::init_singleton()
{
	deserialise_or_init< false >(
		*this,
		fs::path(
			BioEnvironment::singleton().get_factor_synonyms_file()
		),
		boost::bind< void >(
			&EquivalentFactors::init_from_biobase,
			_1
		) 
	);
}



EquivalentFactors::partition_set_t
EquivalentFactors::get_factors_for(BiobaseTablePssmEntry * pssm) const
{
	partition_set_t result;

	const FactorLinkList & factors = pssm->get_factors();

	//for each factor
	for (FactorLinkList::const_iterator f = factors.begin();
		factors.end() != f;
		++f)
	{
		partition_ptr_t partition = find_partition(f->get()->link.entry_idx);
		if (0 != partition)
		{
			result.insert(partition);
		}
	}

	return result;
}



std::string
EquivalentFactors::get_name_for(partition_ptr_t factor)
{
	if (factor->empty())
	{
		throw std::logic_error("Factor partition empty");
	}

	Factor * f = BiobaseDb::singleton().get_entry< FACTOR_DATA >(*(factor->begin()));
	if (0 == f)
	{
		throw std::logic_error("Could not find factor in Biobase");
	}

	return f->get_name();
}







const EquivalentFactorKeyTransformer::result_type &
EquivalentFactorKeyTransformer::operator()( const result_type & partition ) const
{
	return partition;
}



EquivalentFactorKeyTransformer::result_type
EquivalentFactorKeyTransformer::operator()( const TableLink & factor ) const
{
	return EquivalentFactors::singleton().get_partition( factor.entry_idx );
}


BIO_NS_END

