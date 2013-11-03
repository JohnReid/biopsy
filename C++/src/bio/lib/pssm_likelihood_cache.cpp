/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"



//#define BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP

#include "bio/biobase_db.h"
#include "bio/pssm_likelihood_cache.h"
#include "bio/biobase_filter.h"
#include "bio/biobase_data_traits.h"
#include "bio/pssm_likelihood.h"
#include "bio/biobase_match.h"
#include "bio/environment.h"
#include "bio/serialisable.h"

#include <boost/shared_ptr.hpp>
#include <boost/progress.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem/path.hpp>
using namespace boost;
namespace fs = boost::filesystem;

#include <fstream>
using namespace std;




BIO_NS_START

void
PssmLikelihoodCache::populate_from_biobase()
{
	BiobasePssmFilter pssm_filter = BiobasePssmFilter::get_all_pssms_filter();

	//for each site
	for( site_filter_it i = get_sites_begin( pssm_filter );
		get_sites_end( pssm_filter ) != i;
		++i)
	{
		if( is_matchable( i->second ) )
		{
			try
			{
				cout << i->first << ": " << i->second->get_name() << endl;
				get_likelihoods( i->first );
			}
			catch( const std::exception & ex )
			{
				cout << "Could not calculate likelihoods for: " << i->first << ": " << ex.what() << "\n";
			}
			catch( const char * ex )
			{
				cout << "Could not calculate likelihoods for: " << i->first << ": " << ex << "\n";
			}
		}
	}

	//for each matrix
	for( matrix_filter_it i = get_matrices_begin( pssm_filter );
		get_matrices_end( pssm_filter ) != i;
		++i)
	{
		if( is_matchable( i->second ) )
		{
			try
			{
				cout << i->first << ": " << i->second->get_name() << endl;
				get_likelihoods( i->first );
			}
			catch( const std::exception & ex )
			{
				cout << "Could not calculate likelihoods for: " << i->first << ": " << ex.what() << "\n";
			}
			catch( const char * ex )
			{
				cout << "Could not calculate likelihoods for: " << i->first << ": " << ex << "\n";
			}
		}
	}
}


void
PssmLikelihoodCache::init_singleton()
{
	try_to_deserialise< false >(
		*this,
		fs::path(
			BioEnvironment::singleton().get_pssm_likelihoods_cache_file()
		)
	);
}

PssmLikelihoodCache::PssmLikelihoodCache(BiobaseDb & biobase_db)
: biobase_db(biobase_db)
{
}


const BiobaseLikelihoods *
PssmLikelihoodCache::get_likelihoods(const PssmLikelihoodCache::key_t & key)
{
	const BiobaseLikelihoods * result = 0;

	//look for it
	likelihood_map_t::const_iterator l = likelihoods.find(key);
	if (likelihoods.end() == l)
	{
		//we could not find it
		boost::progress_timer timer;
		std::cout << "Calculating pssm likelihoods for: " << key << "\n";

		//calculate the pssm likelihoods
		BiobaseLikelihoods pssm_likelihoods(BioEnvironment::singleton().num_normalisation_quanta);

		pssm_likelihood(
			make_pssm( key ),
			BioEnvironment::singleton().max_pssm_likelihood_map_size,
			pssm_likelihoods );

		//we found it and calculated it so insert it
		std::pair<likelihood_map_t::iterator, bool> insert_result =
			likelihoods.insert( std::make_pair( key, pssm_likelihoods ) );
		assert(insert_result.second); //make sure we did insert it
		result = &(insert_result.first->second);
	}
	else
	{
		//it is already in map
		result = &(l->second); //retrieve it from map
	}

	return result;
}

bool
PssmLikelihoodCache::operator==(const PssmLikelihoodCache & rhs) const
{
	return likelihoods == rhs.likelihoods;
}



BIO_NS_END
