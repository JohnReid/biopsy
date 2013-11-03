/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/xml_builder.h"
#include "bio/x_str.h"

XERCES_CPP_NAMESPACE_USE


BIO_NS_START

using XERCES_CPP_NAMESPACE::DOMDocument;

XmlStartElement::XmlStartElement(const std::string & type, const std::string & ns_uri)
: type(type)
, ns_uri(ns_uri)
{
}

XmlSetTextContent::XmlSetTextContent(const std::string & text)
: text(text)
{
}


XmlBuilder::XmlBuilder(
	DOMDocument * doc,
	DOMElement * current_element )
	: doc( doc )
	, current_element( current_element )
{
}

XmlBuilder &
XmlBuilder::operator<<(const XmlStartElement & x)
{
	DOMElement * child_el = "" != x.ns_uri ? doc->createElementNS(XStr(x.ns_uri), XStr(x.type)) : doc->createElement(XStr(x.type));
	current_element->appendChild(child_el);
	current_element = child_el;
	return *this;
}

XmlBuilder &
XmlBuilder::operator<<(const XmlEndElement & x)
{
	current_element = (DOMElement *)(current_element->getParentNode());

	return *this;
}

XmlBuilder &
XmlBuilder::operator<<(const XmlSetAttribute & x)
{
	if( "" != x.ns_uri ) {
		current_element->setAttributeNS(XStr(x.ns_uri), XStr(x.name), XStr(x.value));
	} else {
		current_element->setAttribute(XStr(x.name), XStr(x.value));
	}
	return *this;
}

XmlBuilder &
XmlBuilder::operator<<(const XmlSetTextContent & x)
{
	current_element->setTextContent(XStr(x.text));
	return *this;
}


XmlBuilder::operator XERCES_CPP_NAMESPACE::DOMElement *()
{
	return current_element;
}


BIO_NS_END
