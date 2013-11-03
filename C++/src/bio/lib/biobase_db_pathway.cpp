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


Pathway::map_t &
BiobaseDb::get_pathways()
{
    return get_deserialise_or_parse< PATHWAY_DATA >( pathways );
}

const Pathway::map_t &
BiobaseDb::get_pathways() const
{
    return get_deserialise_or_parse< PATHWAY_DATA >( pathways );
}





BIO_NS_END
