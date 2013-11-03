/**
@file

Copyright John Reid 2006

*/

#include "biopsy/defs.h"
#include "bio/options.h"


namespace biopsy
{


	void init()
	{
		USING_BIO_NS;

		char * name = const_cast< char * >( "biopsy" );
		BioOptions::singleton().parse( 1, &name );
	}


	std::string get_build()
	{
		const std::string build( 
#ifdef _DEBUG
			"debug"
#else
			"release"
#endif //_DEBUG
			);

		return 
			BIOPSY_MAKE_STRING(
				"0.0.1 biopsy "
				<< build );
	}

}

