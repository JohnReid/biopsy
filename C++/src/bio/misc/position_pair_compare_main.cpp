/* Copyright John Reid 2007
*/

#include "bio-pch.h"



#include "bio/application.h"
#include "bio/matrix_dependencies.h"
USING_BIO_NS

namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;




struct PositionPairCompareApp : Application
{
	PositionPairCompareApp()
	{
		get_options().add_options()
			;
	}

	pssm_source_map_t pssm_source_map;

	int
	task()
	{
		build_all_pssm_sources(pssm_source_map);

		std::copy(pssm_source_map.begin(), pssm_source_map.end(), std::ostream_iterator< pssm_source_map_t::value_type >(cout, "\n"));

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return PositionPairCompareApp().main(argc, argv);
}

