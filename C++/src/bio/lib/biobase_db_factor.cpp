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


Factor::map_t &
BiobaseDb::get_factors()
{
    return get_deserialise_or_parse< FACTOR_DATA >( factors );
}

const Factor::map_t &
BiobaseDb::get_factors() const
{
    return get_deserialise_or_parse< FACTOR_DATA >( factors );
}




BIO_NS_END
