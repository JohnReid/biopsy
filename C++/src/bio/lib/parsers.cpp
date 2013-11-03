/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "DatabaseRefParser.hpp"

BIO_NS_START

bool 
DatabaseRefParser::want_to_parse() const
{
	switch( db )
	{
		//case EMBL_DB:
		case ENSEMBL_DB:
		case ENTREZ_GENE_DB:
		case ENTREZ_PROTEIN_DB:
		case FLYBASE_DB:
		case INPARANOID_DB:
		//case MGI_DB:
		case SWISSPROT_DB:
		case TRANSCOMPEL_DB:
		case TRANSFAC_DB:
		case TRANSPATH_DB:
			return true;
		default:
			return false;
	}
}

BIO_NS_END
