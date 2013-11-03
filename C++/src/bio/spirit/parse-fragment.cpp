#include "bio-pch.h"

#include "parse-table.h"

int
main( int argc, const char * argv[] ) {
	return parse_table< Fragment >( argc, argv );
}
