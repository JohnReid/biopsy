/**
@file

Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"



#include <bio/application.h>
#include <bio/pssm_likelihood_cache.h>
#include <bio/environment.h>
#include <bio/serialisable.h>
USING_BIO_NS;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using namespace std;

/**
 * Estimates the likelihoods of biobase scores (under a binding model) for all PSSMs.
 */
struct CalculatePssmLikelihoodsApp : Application
{
	int task()
	{
		PssmLikelihoodCache::singleton().populate_from_biobase();

		boost::filesystem::path file(
			BioEnvironment::singleton().get_pssm_likelihoods_cache_file()
		);

		std::cout
			<< "Serialising to "
			<< file._BOOST_FS_NATIVE()
			<< "\n"
			;
		serialise< false >(
			PssmLikelihoodCache::singleton(),
			file);

		return 0;
	}
};

int
main(int argc, char * argv [])
{
	return CalculatePssmLikelihoodsApp().main( argc, argv );
}
