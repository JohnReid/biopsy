/* Copyright John Reid 2007
*/

#include "bio-pch.h"


/*
 * Copyright 2002,2004 The Apache Software Foundation.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * $Id: DOMPrint.cpp,v 1.61 2004/09/08 13:55:31 peiyongz Exp $
 */

// ---------------------------------------------------------------------------
//  This sample program invokes the XercesDOMParser to build a DOM tree for
//  the specified input file. It then invokes DOMWriter::writeNode() to
//  serialize the resultant DOM tree back to XML stream.
//
//  Note:
//  Application needs to provide its own implementation of
//		   DOMErrorHandler (in this sample, the DOMPrintErrorHandler),
//		   if it would like to receive notification from the serializer
//		   in the case any error occurs during the serialization.
//
//  Application needs to provide its own implementation of
//		   DOMWriterFilter (in this sample, the DOMPrintFilter),
//		   if it would like to filter out certain part of the DOM
//		   representation, but must be aware that thus may render the
//		   resultant XML stream invalid.
//
//  Application may choose any combination of characters as the
//		   end of line sequence to be used in the resultant XML stream,
//		   but must be aware that thus may render the resultant XML
//		   stream ill formed.
//
//  Application may choose a particular encoding name in which
//		   the output XML stream would be, but must be aware that if
//		   characters, unrepresentable in the encoding specified, appearing
//		   in markups, may force the serializer to terminate serialization
//		   prematurely, and thus no complete serialization would be done.
//
//  Application shall query the serializer first, before set any
//           feature/mode(true, false), or be ready to catch exception if this
//           feature/mode is not supported by the serializer.
//
//  Application needs to clean up the filter, error handler and
//		   format target objects created for the serialization.
//
//   Limitations:
//      1.  The encoding="xxx" clause in the XML header should reflect
//          the system local code page, but does not.
//      2.  Cases where the XML data contains characters that can not
//          be represented in the system local code page are not handled.
//
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//  Includes
// ---------------------------------------------------------------------------
#include "bio/defs.h"

#include "DOMPrintFilter.hpp"
#include "DOMPrintErrorHandler.hpp"


#include <xercesc/util/PlatformUtils.hpp>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMWriter.hpp>
#include <xercesc/dom/DOMWriterFilter.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include <xercesc/util/OutOfMemoryException.hpp>

#include <iostream>
#include <string.h>
#include <stdlib.h>


XERCES_CPP_NAMESPACE_USE





BIO_NS_START

// ---------------------------------------------------------------------------
//
//  dom_print
//
// output_file is null for stdout
//
// ---------------------------------------------------------------------------
void dom_print(DOMNode * doc, const char * output_file)
{
    DOMPrintFilter   *myFilter = 0;

    // get a serializer, an instance of DOMWriter
    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl          = DOMImplementationRegistry::getDOMImplementation(tempStr);
    DOMWriter         *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();

	if (false)
	{
		// even we say to show attribute, but the DOMWriter
		// will not show attribute nodes to the filter as
		// the specs explicitly says that DOMWriter shall
		// NOT show attributes to DOMWriterFilter.
		//
		// so DOMNodeFilter::SHOW_ATTRIBUTE has no effect.
		// same DOMNodeFilter::SHOW_DOCUMENT_TYPE, no effect.
		//
		myFilter = new DOMPrintFilter(DOMNodeFilter::SHOW_ELEMENT   |
										DOMNodeFilter::SHOW_ATTRIBUTE |
										DOMNodeFilter::SHOW_DOCUMENT_TYPE);
		theSerializer->setFilter(myFilter);
	}

	if ( theSerializer->canSetFeature(XMLUni::fgDOMWRTFormatPrettyPrint, true) )
  		theSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, true);
	if (theSerializer->canSetFeature(XMLUni::fgDOMWRTDiscardDefaultContent, true)) {
		theSerializer->setFeature(XMLUni::fgDOMWRTDiscardDefaultContent, true);
	}

	// plug in user's own error handler
    DOMErrorHandler *myErrorHandler = new DOMPrintErrorHandler();
    theSerializer->setErrorHandler(myErrorHandler);

    //
    // Plug in a format target to receive the resultant
    // XML stream from the serializer.
    //
    // StdOutFormatTarget prints the resultant XML stream
    // to stdout once it receives any thing from the serializer.
    //
    XMLFormatTarget *myFormTarget;
    if (output_file)
        myFormTarget = new LocalFileFormatTarget(output_file);
    else
        myFormTarget = new StdOutFormatTarget();

    //
    // do the serialization through DOMWriter::writeNode();
    //
    theSerializer->writeNode(myFormTarget, *doc);

    delete theSerializer;

    //
    // Filter, formatTarget and error handler
    // are NOT owned by the serializer.
    //
    delete myFormTarget;
    delete myErrorHandler;
    delete myFilter;
}

BIO_NS_END
