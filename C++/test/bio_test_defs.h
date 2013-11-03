/**
@file

Copyright John Reid 2013
*/

#ifndef BIO_TEST_DEFS_H_
#define BIO_TEST_DEFS_H_


#ifdef BIO_STANDALONE_TEST
# include <boost/test/included/unit_test.hpp>

/**
 * Define the test suite initialisation function. Used when the test is compiled as a standalone test.
 */
#define BIO_DEFINE_STANDALONE_TEST( name, registration_fn ) \
    test_suite * \
    init_unit_test_suite( int argc, char * argv[] ) \
    { \
        using namespace boost::unit_test; \
        test_suite * ts = BOOST_TEST_SUITE( name ); \
        registration_fn( ts ); \
        return ts; \
    }
#endif //BIO_STANDALONE_TEST


#endif //BIO_TEST_DEFS_H_
