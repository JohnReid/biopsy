/* Copyright John Reid 2007
*/

#include "bio-pch.h"



#include <bio/defs.h>
#include <bio/application.h>
#include <bio/g_soap.h>
USING_BIO_NS

namespace po = boost::program_options;

#include <iostream>
using namespace std;

#include "Generated/soapKEGGBindingProxy.h" // obtain the generated stub
#include "Generated/KEGGBinding.nsmap" // obtain the namespace mapping table


struct KeggApp
	: Application
{
	string request_type;
	bool quiet;
	KEGGBinding kegg_ws; //KEGG web service proxy

	KeggApp()
		: quiet( false )
	{
		get_options().add_options()
		    ("type,t", po::value( &request_type )->default_value( "bget" ), "request type")
		    ("quiet,q", po::bool_switch( &quiet ), "no output")
			;
	}

	void do_query( const std::string & query )
	{
		string result;
		if( ! quiet ) cout << request_type << ": \"" << query << "\"" << endl;
		if( "binfo" == request_type )
		{
			BIO_CHECKED_SOAP_CALL( kegg_ws.ns1__binfo( query, result ) );
		}
		else if( "bget" == request_type )
		{
			BIO_CHECKED_SOAP_CALL( kegg_ws.ns1__bget( query, result ) );
		}
		else if( "bfind" == request_type )
		{
			BIO_CHECKED_SOAP_CALL( kegg_ws.ns1__bfind( query, result ) );
		}
		else if( "btit" == request_type )
		{
			BIO_CHECKED_SOAP_CALL( kegg_ws.ns1__btit( query, result ) );
		}
		else if( "bconv" == request_type )
		{
			BIO_CHECKED_SOAP_CALL( kegg_ws.ns1__bconv( query, result ) );
		}
		else
		{
			throw std::logic_error( BIO_MAKE_STRING( "Unknown request type: " << request_type << ": try one of binfo, bget, bfind, btit, bconv" ) );
		}
		cout << result << endl;
	}

	void check_args()
	{
		if( request_type != "binfo"
			&&
			request_type != "bget"
			&&
			request_type != "bfind"
			&&
			request_type != "btit"
			&&
			request_type != "bconv" )
		{
			throw std::logic_error( BIO_MAKE_STRING( "Unknown request type: " << request_type << ": try one of binfo, bget, bfind, btit, bconv" ) );
		}
			
	}

	int task()
	{
		check_args();

		//initialise runtime environment - once only
		if( ! quiet ) cout << "Initialising SOAP\n";
		g_soap::singleton();

		while( cin )
		{
			if( ! quiet ) cout << "Enter a KEGG " << request_type << " query (or <newline> to exit)" << endl;
			char line[ 2048 ];
			if( ! cin.getline( line, sizeof( line ) - 1 ) || string( "" ) == line )
			{
				break;
			}

			do_query( line );
		}

		return 0;
	}
};


int
main(
	 int argc,
	 char * argv[] )
{
	return KeggApp().main( argc, argv );
}

