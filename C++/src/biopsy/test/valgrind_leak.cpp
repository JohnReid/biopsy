/**
 * Copyright John Reid 2013
 *
 * @file Code to test valgrind return.
 */

#include <iostream>

int
main( int argc, char * argv[] ) {
    int * int_array = new int[ 1 ];
    int_array[0] = 1;
    delete [] int_array;
    std::cout << int_array[0] << "\n"; // <--- this should cause a valgrind "invalid read" error
    return 0;
}
