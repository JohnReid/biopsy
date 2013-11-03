
#ifndef BIO_EXCEPTIONS_H_
#define BIO_EXCEPTIONS_H_

#include "bio/defs.h"

#include <antlr/ANTLRException.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <string>




BIO_NS_START

void translate_dom_exception(const XERCES_CPP_NAMESPACE::DOMException & ex);
void translate_antlr_exception(const antlr::ANTLRException & ex);
void translate_string_exception(const std::string & ex);
void translate_char_exception(const char * ex);

std::string get_last_known_bio_exception();

BIO_NS_END


#endif //BIO_EXCEPTIONS_H_
