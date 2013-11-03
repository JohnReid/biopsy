/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/environment.h"

BIO_NS_START

std::ostream & log_stream()
{
	return *BioEnvironment::singleton().get_log_stream();
}

BIO_NS_END
