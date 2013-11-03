/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/gsl.h"

#include <gsl/gsl_errno.h>

#include <iostream>

BIO_NS_START

void
bio_gsl_error_handler(
	const char * reason,
	const char * file,
	int line,
	int gsl_errno)
{
	using namespace std;

	const std::string msg = BIO_MAKE_STRING("GSL error(" << gsl_errno << "): " << file << "(" << line << "): " << reason);
	cout << msg << endl;
	throw msg;
}

void gsl_init()
{
	gsl_set_error_handler(bio_gsl_error_handler);
}

template <>
RAII<gsl_integration_workspace *>::~RAII()
{
	gsl_integration_workspace_free(resource);
}

BIO_NS_END
