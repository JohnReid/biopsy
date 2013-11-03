
#ifndef BIO_SVG_H_
#define BIO_SVG_H_


#include "bio/defs.h"
#include "bio/x_str.h"

#include <xercesc/dom/DOM.hpp>

#include <string>


BIO_NS_START

typedef std::pair<float_t, float_t> coord_t;

std::string num2str(double num);
//std::string num2str(float_t num);
std::string num2str(size_t num);

void set_attribute(XERCES_CPP_NAMESPACE::DOMElement * el, const XStr & attr, const XStr & value);
void set_attribute(XERCES_CPP_NAMESPACE::DOMElement * el, const XStr & attr, double value);
//void set_attribute(XERCES_CPP_NAMESPACE::DOMElement * el, const XStr & attr, float_t value);
void set_attribute(XERCES_CPP_NAMESPACE::DOMElement * el, coord_t value);

void dom_print(XERCES_CPP_NAMESPACE::DOMNode * doc, const char * output_file = 0);

BIO_NS_END


#endif //BIO_SVG_H_

