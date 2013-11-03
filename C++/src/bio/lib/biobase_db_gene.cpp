/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/biobase_db.h"
#include "bio/biobase_parse_spirit.h"
#include "bio/serialisable.h"
USING_BIO_NS

using namespace boost;
using namespace boost::archive;

using namespace std;



BIO_NS_START



Gene::map_t &
BiobaseDb::get_genes()
{
    return get_deserialise_or_parse< GENE_DATA >( genes );
}

const Gene::map_t &
BiobaseDb::get_genes() const
{
    return get_deserialise_or_parse< GENE_DATA >( genes );
}




BIO_NS_END
