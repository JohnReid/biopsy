
#ifndef BIO_GSOAP_H_
#define BIO_GSOAP_H_

#include "bio/defs.h"
#include "bio/singleton.h"

#include "stdsoap2.h"



#define BIO_CHECKED_SOAP_CALL( x ) if( ( x ) != SOAP_OK ) { throw std::logic_error( g_soap::singleton().get_fault() ); }



BIO_NS_START




/**
Initialises GSOAP library.
*/
struct g_soap
	: Singleton< g_soap >
{
	soap runtime; // gSOAP runtime environment

	~g_soap();

	operator soap & ();
	operator soap * ();

	void init_singleton();

	std::string get_fault();
};


BIO_NS_END

#endif //BIO_GSOAP_H_

