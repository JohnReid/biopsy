
#ifndef BIO_SVG_PSSM_LOGO_H_
#define BIO_SVG_PSSM_LOGO_H_

#include "bio/defs.h"
#include "bio/pssm.h"

#include <xercesc/dom/DOM.hpp>

BIO_NS_START

/** Adds the definitions to the document that are necessary for the pssm logos to work. Call this once per
document */
XERCES_CPP_NAMESPACE::DOMElement *
add_logo_defs(
	XERCES_CPP_NAMESPACE::DOMDocument * doc);


/** Creates an element that visualises a PSSM logo. If a sequence is provided, highlights those bases in the
logo that were in the sequence. */
XERCES_CPP_NAMESPACE::DOMElement *
create_svg_pssm_logo(
	const Pssm & pssm,
	XERCES_CPP_NAMESPACE::DOMDocument * doc,
	const seq_t & seq = "");


BIO_NS_END



#endif //BIO_SVG_PSSM_LOGO_H_
