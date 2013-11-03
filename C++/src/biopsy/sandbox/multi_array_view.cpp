/**
@file

Copyright John Reid 2013

*/



#include <boost/multi_array.hpp>
#include <iostream>


using boost::extents;
using boost::indices;
using namespace std;


typedef boost::multi_array<double, 2> array_t;
typedef array_t::array_view< 2 >::type view;
typedef array_t::const_array_view< 2 >::type const_view;
typedef boost::multi_array_types::index_range range;


const int W = 8;

/// A reverse complement view
view
make_view( array_t a ) {
    return a[ indices[ range( W - 1, -1, -1 ) ] [range( 3, -1, -1 ) ] ];
}

/// A reverse complement view
const_view
make_const_view( array_t a ) {
    return static_cast< const array_t & >( a )[ indices[ range( W - 1, -1, -1 ) ] [range( 3, -1, -1 ) ] ];
}


template< typename Array >
void
print( Array a ) {
    for( array_t::index w = 0; w != W; ++w ) {
        for( array_t::index b = 0; b != 4; ++b ) {
            cout << a[ w ][ b ] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}


int
main( int argc, char * argv[] ) {
    //
    // Initialise our array_t
    //
    double data[] = {
        .7,.1,.1,.1,
        .1,.1,.1,.7,
        .1,.7,.1,.1,
        .1,.1,.7,.1,
        .7,.1,.1,.1,
        .1,.1,.1,.7,
        .1,.7,.1,.1,
        .1,.1,.7,.1,
    };
    array_t a( extents[ W ][ 4 ] );
    a.assign( data, data + W * 4) ;
    print( a );
    print( make_view( a ) );
    print( make_const_view( a ) );

    return 0;
}



