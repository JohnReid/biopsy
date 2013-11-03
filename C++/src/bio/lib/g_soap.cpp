/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/g_soap.h"

BIO_NS_START

g_soap::~g_soap()
{
	soap_destroy( &runtime ); // delete deserialized class instances (for C++ only)
	soap_end( &runtime ); // remove deserialized data and clean up
	soap_done( &runtime ); // detach the gSOAP environment
}

g_soap::operator soap & ()
{
	return runtime;
}

g_soap::operator soap * ()
{
	return &runtime;
}

void
g_soap::init_singleton()
{
	if( soap_check_state( &runtime ) )
	{
		throw std::logic_error( "GSOAP already initialised" );
	}

	soap_init( &runtime ); // initialize runtime environment (only once)
}

std::string
g_soap::get_fault()
{
	const char *c, *v = NULL, *s, **d;
	d = soap_faultcode( &runtime );
	if( ! *d )
	{
		soap_set_fault( &runtime );
	}
	c = *d;
	if( runtime.version == 2 )
	{
		v = *soap_faultsubcode( &runtime );
	}
	s = *soap_faultstring( &runtime );
	d = soap_faultdetail( &runtime );

	std::stringstream str;
	return
		BIO_MAKE_STRING(
			"SOAP 1."
			<< int( runtime.version )
			<< ": "
			<< c
			<< " ["
			<< ( v ? v : "no subcode" )
			<< "] \""
			<< ( s ? s : "[no reason]" )
			<< "\" "
			<< ( d && *d ? *d : "[no detail]" ) );
}


BIO_NS_END
