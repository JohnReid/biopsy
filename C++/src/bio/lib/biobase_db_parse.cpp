/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
USING_BIO_NS

using namespace boost;
using namespace boost::archive;

using namespace std;



BIO_NS_START




BiobaseTablePssmEntry *
BiobaseDb::get_pssm_entry( const TableLink & link ) const
{
	switch(link.table_id)
	{
		case SITE_DATA:
			//site
			return get_entry< SITE_DATA >( link ) ;

		case MATRIX_DATA:
			//matrix
			return get_entry< MATRIX_DATA >( link );
		default: throw std::logic_error( "Bad table index" );
	}
}

BiobaseTableEntry *
BiobaseDb::get_entry(const TableLink & link) const
{
	switch(link.table_id) {
		case MOLECULE_DATA:
			//molecule
			return get_entry<MOLECULE_DATA>(link);

		case PATHWAY_DATA:
			//pathway
			return get_entry<PATHWAY_DATA>(link);

		case FRAGMENT_DATA:
			//fragment
			return get_entry<FRAGMENT_DATA>(link);

		case GENE_DATA:
			//gene
			return get_entry<GENE_DATA>(link);

		case FACTOR_DATA:
			//factor
			return get_entry<FACTOR_DATA>(link);

		case SITE_DATA:
			//site
			return get_entry<SITE_DATA>(link);

		case MATRIX_DATA:
			//matrix
			return get_entry<MATRIX_DATA>(link);
		default: throw std::logic_error( "Bad table index" );
	}
}



void
BiobaseDb::load_all() const 
{
	get_matrices();
	get_sites();
	get_factors();
	get_fragments();
	get_genes();
	get_compels();
	get_evidences();
//	get_pathways();
//	get_molecules();
}




BIO_NS_END
