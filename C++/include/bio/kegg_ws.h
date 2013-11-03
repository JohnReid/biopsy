

#ifndef BIO_KEGG_WS_H_
#define BIO_KEGG_WS_H_

#include "bio/singleton.h"
#include "bio/bimap.h"

#include "Generated/soapKEGGBindingProxy.h" // obtain the generated stub
#include "Generated/KEGGBinding.nsmap" // obtain the namespace mapping table

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>

BIO_NS_START



template< 
	typename VertexKey = std::string,
	typename EdgeKey = std::string
>
struct gene_network
{
	typedef
		boost::adjacency_list <
			boost::vecS,
			boost::vecS,
			boost::undirectedS,
			boost::property <
				boost::vertex_name_t,
				VertexKey >,
			boost::property <
				boost::edge_name_t,
				EdgeKey >
		>
		graph;

	typedef
		typename boost::property_map <
			graph,
			boost::vertex_name_t >::type
		id_map;

	typedef
		typename boost::property_map <
			graph,
			boost::edge_name_t >::type
		pathway_map;

	typedef
		typename boost::graph_traits < graph >::vertex_descriptor
		vertex;

	typedef
		typename boost::graph_traits < graph >::edge_descriptor
		edge;

	typedef
		typename bimap<
			std::string,
			vertex >::unique
		gene_id_vertex_map;

	graph g;
	gene_id_vertex_map ids;
	id_map id;
	pathway_map pathway;

	gene_network()
	{
		id = boost::get( boost::vertex_name, g );
		pathway = boost::get( boost::edge_name, g );
	}

	vertex get_or_insert_vertex( const VertexKey & key )
	{
		typename gene_id_vertex_map::iterator pos = ids.find( key );
		if( ids.end() == pos )
		{
			pos = ids.insert( std::make_pair( key, add_vertex( g ) ) ).first;
			id[ pos->second ] = key;
		}

		return pos->second;
	}

	edge get_or_insert_edge( vertex u, vertex v, const EdgeKey & key )
	{
		edge e;
		bool inserted;
		boost::tie( e, inserted ) = add_edge( u, v, g );
		if( inserted )
		{
			pathway[ e ] = key;
		}

		return e;
	}


};


/**
Matrix should conform to the boost.uBLAS matrix concept.

Network should be a gene network.
*/
template<
	typename Graph,
	typename Matrix >
struct build_adjacency_matrix_graph_visitor
{
	typedef Graph graph;
	typedef typename graph::vertex_descriptor vertex;
	typedef typename graph::edge_descriptor edge;
	typedef Matrix matrix;
	typedef bimap< vertex, unsigned > idx_map;

	unsigned idx;
	matrix & m;
	typename idx_map::unique indices;

	build_adjacency_matrix_graph_visitor( matrix & m )
		: m( m )
		, idx( 0 )
	{
	}

	/** define the mapping between the vertices and their indices */
	void initialize_vertex( const vertex & v, const graph & g )
	{
		indices.insert( idx_map::value_type( v, idx++ ) );
	}

	/** Put the edge in the adjacency matrix. */
	void examine_edge( const edge & e, const graph & g )
	{
		if( m.size1() != idx || m.size2() != idx )
		{
			m.resize( idx, idx, false );
			m.clear();
		}

		vertex u, v;
		boost::tie( u, v ) = boost::incident( e, g );
		const unsigned i = indices.find( u )->second;
		const unsigned j = indices.find( v )->second;
		if( 0 == m( i, j ) )
		{
			m.insert_element( i, j, 1 );
		}
		else
		{
			int i = 0; //already in the matrix
		}
	}

	void start_vertex( const vertex & v, const graph & g ) { }
	void discover_vertex( const vertex & v, const graph & g ) { }
	void examine_vertex( const vertex & v, const graph & g ) { }	
	void tree_edge( const edge & e, const graph & g ) { }
	void back_edge( const edge & e, const graph & g ) { }
	void forward_or_cross_edge( const edge & e, const graph & g ) { }
	void finish_vertex( const vertex & v, const graph & g ) { }
};


/**
Matrix should conform to the boost.uBLAS matrix concept.

Network should be a gene network.
*/
template<
	typename Graph,
	typename Matrix >
void
fill_adjacency_matrix(
	const Graph & g,
	Matrix & a )
{
	using namespace boost;
	using namespace boost::graph;

#if 0 //using visitor
	depth_first_search(
		g,
		visitor( build_adjacency_matrix_graph_visitor< Graph, Matrix >( a ) ) );
#endif


	//associate an index from 0,...,n with each vertex
	typedef bimap< typename graph_traits< Graph >::vertex_descriptor, unsigned > idx_map;
	typename idx_map::unique vertex_indices;
	unsigned idx = 0;
	typename graph_traits< Graph >::vertex_iterator v_begin, v_end;
	tie( v_begin, v_end ) = vertices( g );
	for( ; v_begin != v_end; ++v_begin, ++idx )
	{
		vertex_indices.insert( typename idx_map::value_type( *v_begin, idx ) );
	}

	//resize and clear matrix
	a.resize( idx, idx, false );
	a.clear();

	//visit each edge to insert into matrix
	typename graph_traits< Graph >::edge_iterator e_begin, e_end;
	tie( e_begin, e_end ) = edges( g );
	for( ; e_begin != e_end; ++e_begin )
	{
		typename graph_traits< Graph >::vertex_descriptor u, v;
		boost::tie( u, v ) = boost::incident( *e_begin, g );
		const unsigned i = vertex_indices.find( u )->second;
		const unsigned j = vertex_indices.find( v )->second;
		a( i, j ) = 1;
		a( j, i ) = 1;
	}
}



/**
The KEGG web service.
*/
struct kegg_ws
	: Singleton< kegg_ws >
{
	KEGGBinding proxy; //KEGG web service proxy

	void init_singleton()
	{
		//initialise runtime environment - once only
		g_soap::singleton();
	}
};


template< typename OutputIt >
void
insert_pathways_from_organism(
	const std::string & organism,
	OutputIt output_it )
{
	//get the pathways from the web service
	ns1__list_USCOREpathwaysResponse results;
	BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__list_USCOREpathways( organism, results ) );
	for( int i = 0; i < results.return_->__size; ++i )
	{
		*output_it =
			std::make_pair(
			//typename OutputIt::value_type(
				( *( results.return_->__ptr + i ) )->entry_USCOREid,
				( *( results.return_->__ptr + i ) )->definition );
		++output_it;
	}
}



template< typename OutputIt >
void
insert_databases(
	OutputIt output_it )
{
	//get the databases from the web service
	ns1__list_USCOREdatabasesResponse results;
	BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__list_USCOREdatabases( results ) );
	for( int i = 0; i < results.return_->__size; ++i )
	{
		*output_it = 
			std::make_pair(
				( *( results.return_->__ptr + i ) )->entry_USCOREid,
				( *( results.return_->__ptr + i ) )->definition );
		++output_it;
	}
}




template< typename OutputIt >
void
insert_organisms(
	OutputIt output_it )
{
	//get the organisms from the web service
	ns1__list_USCOREorganismsResponse results;
	BIO_CHECKED_SOAP_CALL( kegg_ws::singleton().proxy.ns1__list_USCOREorganisms( results ) );
	for( int i = 0; i < results.return_->__size; ++i )
	{
		*output_it = 
			std::make_pair(
				( *( results.return_->__ptr + i ) )->entry_USCOREid,
				( *( results.return_->__ptr + i ) )->definition );
		++output_it;
	}
}

template< typename OutputIt >
void
insert_genes_from_pathway(
	const std::string & pathway,
	OutputIt output_it )
{
	//get the databases from the web service
	ns1__get_USCOREgenes_USCOREby_USCOREpathwayResponse results;
	BIO_CHECKED_SOAP_CALL(
		kegg_ws::singleton().proxy.ns1__get_USCOREgenes_USCOREby_USCOREpathway(
			pathway,
			results ) );
	for( int i = 0; i < results.return_->__size; ++i )
	{
		*output_it = **( results.return_->__ptr + i );
		++output_it;
	}
}

unsigned
get_num_genes_in_organism( const std::string & organism )
{
	int num_genes;
	BIO_CHECKED_SOAP_CALL(
		kegg_ws::singleton().proxy.ns1__get_USCOREnumber_USCOREof_USCOREgenes_USCOREby_USCOREorganism(
			organism,
			num_genes ) );
	return num_genes;
}

template< typename GeneInsIt >
void
insert_genes_from_organism(
	const std::string & organism,
	GeneInsIt gene_ins_it )
{
	int num_genes;
	BIO_CHECKED_SOAP_CALL(
		kegg_ws::singleton().proxy.ns1__get_USCOREnumber_USCOREof_USCOREgenes_USCOREby_USCOREorganism(
			organism,
			num_genes ) );

	ns1__get_USCOREgenes_USCOREby_USCOREorganismResponse results;
	BIO_CHECKED_SOAP_CALL(
		kegg_ws::singleton().proxy.ns1__get_USCOREgenes_USCOREby_USCOREorganism(
			organism,
			1,
			num_genes,
			results ) );

	if( results.return_->__size != num_genes )
	{
		throw
			std::logic_error(
				BIO_MAKE_STRING(
					"Unexpected # genes; "
					<< results.return_->__size
					<< " instead of "
					<< num_genes ) );
	}

	for( int i = 0; results.return_->__size != i; ++i )
	{
		*gene_ins_it = **( results.return_->__ptr + i );
		++gene_ins_it;
	}
}



unsigned get_num_pathways_for( const std::string & gene_id_1, const std::string & gene_id_2 )
{
	//build input
	ArrayOfstring genes;
	genes.__size = 2;
	genes.__ptr = new std::string *[ 2 ];
	genes.__ptr[ 0 ] = &const_cast< std::string & >( gene_id_1 );
	genes.__ptr[ 1 ] = &const_cast< std::string & >( gene_id_2 );

	ns1__get_USCOREpathways_USCOREby_USCOREgenesResponse results;
	BIO_CHECKED_SOAP_CALL(
		kegg_ws::singleton().proxy.ns1__get_USCOREpathways_USCOREby_USCOREgenes(
			&genes,
			results ) );

	return results.return_->__size;
}



void
build_network_for_organism(
	const std::string & organism,
	gene_network<> & network,
	bool verbose = false )
{
	typedef std::map< std::string, std::string > map;
	map pathways;
	insert_pathways_from_organism( organism, std::inserter( pathways, pathways.begin() ) );

	if( verbose )
	{
		std::cout
			<< "Got "
			<< pathways.size()
			<< " pathways for organism \""
			<< organism
			<< "\"\n";

		std::cout
			<< "\""
			<< organism
			<< "\" has "
			<< get_num_genes_in_organism( organism )
			<< " genes\n";
	}

	BOOST_FOREACH( const map::value_type & pathway, pathways )
	{
		std::vector< std::string > genes;
		insert_genes_from_pathway( pathway.first, std::back_inserter( genes ) );

		for( unsigned i1 = 0; genes.size() != i1; ++i1 )
		{
			gene_network<>::vertex v1 = network.get_or_insert_vertex( genes[ i1 ] );

			for( unsigned i2 = i1 + 1; genes.size() != i2; ++i2 )
			{
				gene_network<>::vertex v2 = network.get_or_insert_vertex( genes[ i2 ] );

				gene_network<>::edge u = network.get_or_insert_edge( v1, v2, pathway.first );
			}
		}
	}

	if( verbose )
	{
		const unsigned num_vertices = boost::num_vertices( network.g );
		const unsigned num_edges = boost::num_edges( network.g );

		std::cout
			<< "Network has "
			<< num_vertices
			<< " vertices and "
			<< num_edges
			<< " edges: "
			<< 200 * num_edges / num_vertices * num_vertices
			<< "% sparsity\n";
	}
}






BIO_NS_END

#endif //BIO_KEGG_WS_H_
