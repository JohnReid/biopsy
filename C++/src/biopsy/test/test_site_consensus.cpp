/**
 * Copyright John Reid 2010, 2013
 *
 * @file Iterates through all the sites in TRANSFAC and gets their IUPAC codes.
 */

#define BOOST_TEST_MODULE site_consensus
#include <boost/test/unit_test.hpp>

#include <biopsy/init.h>
#include <bio/biobase_db.h>

USING_BIO_NS;

BOOST_AUTO_TEST_CASE( test_site_consensus )
{
    using namespace biopsy;

    init();

    BOOST_FOREACH( Site::map_t::value_type const & site, BiobaseDb::singleton().get_sites() ) {
        std::cout
            << site.first << ": "
            << site.second->sequence
            << "\n"
            ;
        // only look at sequences without 'U's
        if( std::string::npos == site.second->sequence.find_first_of( "uUX." ) ) {
            std::cout
                << "     -> "
                ;
            BOOST_FOREACH( char s, site.second->sequence ) {
                try {
                    std::cout << IupacCode( s );
                } catch( const std::logic_error & e ) {
                    std::cout << "X";
                }
            }
            std::cout << "\n";
        }
    }
}
