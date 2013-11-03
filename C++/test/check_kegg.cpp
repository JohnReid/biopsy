

#include <bio/defs.h>
#include <bio/g_soap.h>
#include <bio/kegg_ws.h>
USING_BIO_NS

#include <boost/foreach.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/progress.hpp>
#undef max
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;



#define VERBOSE_CHECKING


static const std::string organism = "bpn";

void
check_organisms()
{
	cout << "******* check_organisms()" << endl;

	typedef std::map< std::string, std::string > map;
	map organisms;
	insert_organisms( std::inserter( organisms, organisms.begin() ) );

#ifdef VERBOSE_CHECKING
	std::cout << "Organism - genes - pathways - desc:\n";
	BOOST_FOREACH( const map::value_type & value, organisms )
	{
		map pathways;
		insert_pathways_from_organism( value.first, std::inserter( pathways, pathways.begin() ) );

		std::cout
			<< value.first
			<< " - " 
			<< get_num_genes_in_organism( value.first )
			<< " - " 
			<< pathways.size()
			<< " - " 
			<< value.second
			<< "\n";
	}

#endif //VERBOSE_CHECKING
}



void
check_get_genes_by_organism()
{
	cout << "******* check_get_genes_by_organism()" << endl;

	std::vector< std::string > genes;
	insert_genes_from_organism( organism, back_inserter( genes ) );
}



void
check_build_gene_network()
{
	cout << "******* check_build_gene_network()" << endl;

	{
		gene_network<> network;
		build_network_for_organism( organism, network, true );

		using namespace boost::numeric::ublas;
		compressed_matrix< unsigned > A;
		fill_adjacency_matrix( network.g, A );
	}
}



void
check_compressed_matrix()
{
	cout << "******* check_gene_adjacency()" << endl;

    using namespace boost::numeric::ublas;

	compressed_matrix< unsigned > test( 10, 10 );
	test( 0, 0 ) = test( 1, 1 ) = 1;
	test( 2, 2 ) = test( 3, 3 ) = 0;
}



void
check_kegg_wsdl()
{
	cout << "******* check_kegg_wsdl()" << endl;

	typedef std::map< std::string, std::string > map;

	{
		map databases;
		insert_databases( std::inserter( databases, databases.begin() ) );

#ifdef VERBOSE_CHECKING
		std::cout << "Databases:\n";
		BOOST_FOREACH( const map::value_type & value, databases )
		{
			std::cout << value.first << " - " << value.second << "\n";
		}
		std::cout << "\n";
#endif //VERBOSE_CHECKING
	}

	{
		map organisms;
		insert_organisms( std::inserter( organisms, organisms.begin() ) );

#ifdef VERBOSE_CHECKING
		std::cout << "Organisms:\n";
		BOOST_FOREACH( const map::value_type & value, organisms )
		{
			std::cout << value.first << " - " << value.second << "\n";
		}
		std::cout << "\n";
#endif //VERBOSE_CHECKING
	}

	{
		map pathways;
		insert_pathways_from_organism( organism, std::inserter( pathways, pathways.begin() ) );

#ifdef VERBOSE_CHECKING
		std::cout << "Pathways for \"" << organism << "\" :\n";
		BOOST_FOREACH( const map::value_type & value, pathways )
		{
			std::cout << value.first << " - " << value.second << "\n";
		}
		std::cout << "\n";
#endif //VERBOSE_CHECKING
	}

	{
		string result;
		cout << "bget: \"path:ko04330\"" << endl;
		BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__bget( "path:ko04330", result ) );
		cout << "result: " << result << endl;
	}

	{
		string result;
		cout << "bfind: \"path:ko04330\"" << endl;
		BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__bfind( "path:ko04330", result ) );
		cout << "result: " << result << endl;
	}

	{
		cout << "get_genes_by_pathway: \"path:ko04330\"" << endl;
		ns1__get_USCOREgenes_USCOREby_USCOREpathwayResponse result;
		BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__get_USCOREgenes_USCOREby_USCOREpathway( "path:eco00020", result ) );
		cout << "result: " << endl;
	}

#ifdef VERBOSE_CHECKING
	while( cin )
	{
		cout << "Enter a KEGG bget query or <newline> to stop" << endl;
		char line[ 2048 ];
		if( ! cin.getline( line, sizeof( line ) - 1 ) || string( "" ) == line )
		{
			break;
		}

		string result;
		cout << "bget: \"" << line << "\"" << endl;
		BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__bget( line, result ) );
		cout << "result: " << result << endl;

	}
#endif //VERBOSE_CHECKING

#ifdef VERBOSE_CHECKING
	while( cin )
	{
		cout << "Enter a KEGG bfind query or <newline> to stop" << endl;
		char line[ 2048 ];
		if( ! cin.getline( line, sizeof( line ) - 1 ) || string( "" ) == line )
		{
			break;
		}

		string result;
		cout << "bfind: \"" << line << "\"" << endl;
		BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__bfind( line, result ) );
		cout << "result: " << result << endl;

	}
#endif //VERBOSE_CHECKING
}


void
register_kegg_tests( boost::unit_test::test_suite * test )
{
	//test->add( BOOST_TEST_CASE( &check_organisms ), 0);
	test->add( BOOST_TEST_CASE( &check_compressed_matrix ), 0);
	test->add( BOOST_TEST_CASE( &check_build_gene_network ), 0);
	//test->add( BOOST_TEST_CASE( &check_kegg_wsdl ), 0);
	//test->add( BOOST_TEST_CASE( &check_get_genes_by_organism ), 0);
}
