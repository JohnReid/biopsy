

#include "bio/defs.h"
#include "bio/environment.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;

#include <wsdlparser/WsdlInvoker.h>
using namespace WsdlPull;


//#define VERBOSE_CHECKING




void
check_wsdl()
{
	cout << "******* check_wsdl()" << endl;

	WsdlInvoker invoker;
	if( ! invoker.setWSDLUri( BioEnvironment::singleton().get_kegg_wsdl_uri() ) )
	{
		throw std::logic_error( invoker.errors() );
	}

#ifdef VERBOSE_CHECKING
	invoker.setVerbose(true);   
#endif

	if( ! invoker.setOperation( "bgetRequest" ) )
	{
		throw std::logic_error( "Error calling bgetRequest" );
	}

	if( ! invoker.status() )
	{
		throw std::logic_error( "Bad invoker status" );
	}

	Schema::Type t;
	const std::string pathway_id = "ko04350";
	if( ! invoker.setValue( "QuoteTicker", ( void* )( &pathway_id ) ) )
	{
		throw std::logic_error( BIO_MAKE_STRING( "Incorrect input value: " << pathway_id ) );
	}

	if (! invoker.invoke() )
	{
		throw std::logic_error( "Couldn't invoke the web service" );
	}
	else
	{
		void * val = invoker.getValue( "bgetResponse", t );
		std::string * string_ptr = static_cast< std::string * >( val );

		std::cout
			<< pathway_id
			<< "-->"
			<< *string_ptr
			<< std::endl;
	}
}



void
register_wsdl_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_wsdl ), 0);
}
