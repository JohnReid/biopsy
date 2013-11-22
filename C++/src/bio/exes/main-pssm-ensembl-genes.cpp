/* Copyright John Reid 2006, 2007
*/

#include "bio-pch.h"

#include "bio/application.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/biobase_filter.h"

#include <boost/foreach.hpp>
#include <boost/range.hpp>

USING_BIO_NS;
using namespace boost;
using namespace std;


/**
Prints a list of PSSMS with the ensembl genes associated with them in csv format. 
*/
struct AliasesApp : Application
{
	unsigned num_not_found_factors;
	unsigned num_not_found_genes;
	unsigned num_no_gene_link;
	unsigned num_pssms_examined;
	unsigned num_genes_found;

	AliasesApp()
		: num_not_found_factors( 0 )
		, num_not_found_genes( 0 )
		, num_no_gene_link( 0 )
		, num_pssms_examined( 0 )
		, num_genes_found( 0 )
	{
		get_options().add_options()
			//( "input,i", po::value( &input )->default_value( "" ), "serialised statistics" )
			;
	}

	template< typename Range >
	void aliases( const Range & range )
	{
		BOOST_FOREACH( const typename Range::value_type & entry, range )
		{
			++num_pssms_examined;

			cout << entry.first;
			BOOST_FOREACH( FactorLinkPtr factor_link, entry.second->get_factors() )
			{
				if( FACTOR_DATA != factor_link->link.table_id )
				{
					throw std::logic_error( 
						BIO_MAKE_STRING( 
							"Factor link does not point to factor data: "
							<< factor_link->link.table_id
							<< " for "
							<< entry.first ) );
				}

				Factor * factor = BiobaseDb::singleton().get_entry< FACTOR_DATA >( factor_link->link );
				if( 0 == factor )
				{
					++num_not_found_factors;
					continue;
				}
				
				if( UNKNOWN_DATA == factor->gene.table_id )
				{
					++num_no_gene_link;
					continue;
				}

				if( GENE_DATA != factor->gene.table_id )
				{
					throw std::logic_error( 
						BIO_MAKE_STRING( 
							"Gene link does not point to gene data: "
							<< factor->gene.table_id
							<< " for "
							<< entry.first ) );
				}

				Gene * gene = BiobaseDb::singleton().get_entry< GENE_DATA >( factor->gene );
				if( 0 == gene )
				{
					++num_not_found_genes;
					continue;
				}

				set< string > accessions;
				BOOST_FOREACH( const DatabaseRefVec::value_type & ref, gene->database_refs )
				{
					if( ENSEMBL_DB == ref.db )
					{
						accessions.insert( BIO_MAKE_STRING ( ref ) );
					}
				}

				BOOST_FOREACH( const string & acc, accessions )
				{
					++num_genes_found;
					cout << "," << acc;
				}
			}
			cout << "\n";
		}
	}

	int task()
	{
		//make sure all is loaded
		BiobaseDb::singleton().get_matrices();
		BiobaseDb::singleton().get_sites();
		BiobaseDb::singleton().get_factors();
		BiobaseDb::singleton().get_genes();
	
		BiobasePssmFilter filter = BiobasePssmFilter::get_all_pssms_filter();
		aliases( 
			make_iterator_range( 
				get_matrices_begin( filter ), 
				get_matrices_end( filter ) ) );
	
		aliases( 
			make_iterator_range( 
				get_sites_begin( filter ), 
				get_sites_end( filter ) ) );

		cout
			<< num_pssms_examined
			<< " pssms examined\n"
			<< num_not_found_factors
			<< " factors not found\n"
			<< num_no_gene_link
			<< " factors without a gene link\n"
			<< num_genes_found
			<< " genes found\n"
			<< num_not_found_genes
			<< " genes not found\n";


		return 0;
	}
};




int
main(
	int argc,
	char * argv[] )
{
	return AliasesApp().main( argc, argv );
}



