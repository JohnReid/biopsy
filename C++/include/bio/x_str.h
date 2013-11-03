#ifndef BIO_X_STR_H_
#define BIO_X_STR_H_

#include "bio/defs.h"

#include <xercesc/util/XMLString.hpp>

#include <iostream>
#include <string>



BIO_NS_START

class XStr
{
public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    XStr(const char* const toTranscode)
    {
        // Call the private transcoding method
        fUnicodeForm = XERCES_CPP_NAMESPACE::XMLString::transcode(toTranscode);
    }

	XStr(const std::string & toTranscode)
    {
        // Call the private transcoding method
        fUnicodeForm = XERCES_CPP_NAMESPACE::XMLString::transcode(toTranscode.c_str());
    }

    ~XStr()
    {
        XERCES_CPP_NAMESPACE::XMLString::release(&fUnicodeForm);
    }


    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const XMLCh* unicodeForm() const
    {
        return fUnicodeForm;
    }

	operator const XMLCh * () const { return fUnicodeForm; }

private :
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fUnicodeForm
    //      This is the Unicode XMLCh format of the string.
    // -----------------------------------------------------------------------
    XMLCh*   fUnicodeForm;
};





// ---------------------------------------------------------------------------
//  This is a simple class that lets us do easy (though not terribly efficient)
//  trancoding of XMLCh data to local code page for display.
// ---------------------------------------------------------------------------
class StrX
{
public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    StrX(const XMLCh* const toTranscode)
    {
        // Call the private transcoding method
        fLocalForm = XERCES_CPP_NAMESPACE::XMLString::transcode(toTranscode);
    }

    ~StrX()
    {
        XERCES_CPP_NAMESPACE::XMLString::release(&fLocalForm);
    }


    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const char* localForm() const
    {
        return fLocalForm;
    }

	operator const char * () const {
		return fLocalForm;
	}

private :
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fLocalForm
    //      This is the local code page form of the string.
    // -----------------------------------------------------------------------
    char*   fLocalForm;
};


inline XERCES_STD_QUALIFIER ostream& operator<<(XERCES_STD_QUALIFIER ostream& target, const StrX& toDump)
{
    target << toDump.localForm();
    return target;
}

BIO_NS_END


#endif //BIO_X_STR_H_

