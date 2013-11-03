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



Matrix::map_t &
BiobaseDb::get_matrices()
{
    return get_deserialise_or_parse< MATRIX_DATA >( matrices );
}

const Matrix::map_t &
BiobaseDb::get_matrices() const
{
    return get_deserialise_or_parse< MATRIX_DATA >( matrices );
}




BIO_NS_END
