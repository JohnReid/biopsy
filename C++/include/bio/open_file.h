#ifndef BIO_OPEN_FILE_H_
#define BIO_OPEN_FILE_H_

#include "bio/defs.h"

#include <boost/filesystem/path.hpp>

#ifdef _MSC_VER
# pragma comment( lib, "shell32" )
#endif

BIO_NS_START

void
open_file(const boost::filesystem::path & file);

BIO_NS_END

#endif //BIO_OPEN_FILE_H_
