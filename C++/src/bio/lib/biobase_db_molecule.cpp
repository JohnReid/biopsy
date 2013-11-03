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



Molecule::map_t &
BiobaseDb::get_molecules()
{
    return get_deserialise_or_parse< MOLECULE_DATA >( molecules );
}

const Molecule::map_t &
BiobaseDb::get_molecules() const
{
    return get_deserialise_or_parse< MOLECULE_DATA >( molecules );
}




BIO_NS_END
