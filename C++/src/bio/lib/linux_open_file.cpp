/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/open_file.h"

#include <boost/filesystem/path.hpp>

BIO_NS_START

void
open_file(const boost::filesystem::path & file)
{
	throw std::logic_error( "open_file() not implemented on linux" );
}

BIO_NS_END

