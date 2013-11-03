
#ifndef BIO_XML_BUILDER_H_
#define BIO_XML_BUILDER_H_

#include "bio/defs.h"
#include "bio/x_str.h"

#include <xercesc/dom/DOM.hpp>

#include <string>


BIO_NS_START


struct XmlStartElement
{
	XmlStartElement(const std::string & type, const std::string & ns_uri = "");
	std::string type;
	std::string ns_uri;
};

struct XmlEndElement
{
};

struct XmlSetAttribute
{
	template <class Value>
	XmlSetAttribute(const std::string & name, const Value & value, const std::string & ns_uri = "");

	std::string name;
	std::string value;
	std::string ns_uri;
};

struct XmlSetTextContent
{
	XmlSetTextContent(const std::string & text);

	std::string text;
};

class XmlBuilder
{
public:
	XmlBuilder(
		XERCES_CPP_NAMESPACE::DOMDocument * doc,
		XERCES_CPP_NAMESPACE::DOMElement * current_element );

	XmlBuilder & operator<<(const XmlStartElement &);
	XmlBuilder & operator<<(const XmlEndElement &);
	XmlBuilder & operator<<(const XmlSetAttribute &);
	XmlBuilder & operator<<(const XmlSetTextContent &);

	operator XERCES_CPP_NAMESPACE::DOMElement *();

	XERCES_CPP_NAMESPACE::DOMDocument * doc;
	XERCES_CPP_NAMESPACE::DOMElement * current_element;
};

/*
template <>
inline
XmlSetAttribute::XmlSetAttribute(const std::string & name, const std::string & value, const std::string & ns_uri)
	: name(name)
	, value(value)
	, ns_uri(ns_uri)
{
}
*/

template <class Value>
XmlSetAttribute::XmlSetAttribute(const std::string & name, const Value & value, const std::string & ns_uri)
	: name(name)
	, value(BIO_MAKE_STRING(value))
	, ns_uri(ns_uri)
{
}

	
BIO_NS_END



#endif //BIO_XML_BUILDER_H_
