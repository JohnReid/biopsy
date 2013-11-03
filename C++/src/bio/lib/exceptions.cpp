/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/exceptions.h"


using namespace antlr;

XERCES_CPP_NAMESPACE_USE

#include <iostream>
using namespace std;


BIO_NS_START

static string last_known_exception;

std::string get_last_known_bio_exception()
{
	std::string result = last_known_exception;
	last_known_exception.clear();
	return result;
}

void translate_dom_exception(const DOMException & ex)
{
	last_known_exception = BIO_MAKE_STRING(ex.msg);
	cerr << last_known_exception << endl;
}

void translate_antlr_exception(const ANTLRException & ex)
{
	last_known_exception = ex.toString();
	cerr << last_known_exception << endl;
}

void translate_string_exception(const string & ex)
{
	last_known_exception = ex;
	cerr << last_known_exception << endl;
}

void translate_char_exception(const char * ex)
{
	last_known_exception = ex;
	cerr << last_known_exception << endl;
}

BIO_NS_END
