#ifndef BIO_BIOBASE_PARSE_SPIRIT_H_
#define BIO_BIOBASE_PARSE_SPIRIT_H_

#include "bio/defs.h"

//
// Turn off uninitialised warnings in boost code.
//
#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wuninitialized"
#endif //__GNUC__

#include "bio/biobase_data_traits.h"
#include "bio/biobase_db.h"
#include "bio/serialisable.h"
#include <bio/spirit/transfac_qi.h>

#include <boost/spirit/include/support_istream_iterator.hpp>
//#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

//
// Turn warnings back on
//
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif //__GNUC__


#include <fstream>

BIO_NS_START
namespace spirit {


template< TransData type >
void
parse(
    typename DataTraits< type >::entry_t::map_t & map,
    const std::string & biobase_file = DataTraits< type >::get_biobase_file()
)
{
    using namespace std;
    using namespace BIO_NS::spirit;

    typedef typename DataTraits< type >::entry_t entry_t;
    typedef typename entry_t::map_t map_t;

    ifstream input( biobase_file.c_str() );
    if( ! input ) {
        throw std::logic_error(
            BIO_MAKE_STRING( "Could not open file: \"" << biobase_file << "\"" ) );
    }
    input.unsetf( std::ios::skipws ); // disable skipping of whitespace

    // iterate over stream input using boost::spirit forward iterator for streams
    typedef boost::spirit::istream_iterator base_iterator_type;
    base_iterator_type input_begin( input );

    // wrap forward iterator with position iterator, to record the position
    typedef boost::spirit::line_pos_iterator< base_iterator_type > iterator_type;

    iterator_type begin( input_begin );
    iterator_type end;

    try
    {
        typedef typename get_entry_parser< entry_t, iterator_type >::type parser;
        biobasetable_parser< iterator_type, parser > table;

        using namespace qi;
        using qi::ascii::blank;

        boost::timer timer;
        iterator_type i = begin;
        qi::parse(
            i,
            end,
            table,
            map
        );
        std::cout
            << "Parsed " << map.size()
            << " entries from \""
            << biobase_file << "\" - "
            << timer.elapsed() << "s\n";

        if( end != i ) {
            throw std::logic_error(
                BIO_MAKE_STRING( "Did not parse all input: Got to line "
                    << boost::spirit::get_line( i )
                    << ", column " << boost::spirit::get_column( begin, i ) ) );
        }

    } catch ( qi::expectation_failure< iterator_type > const & e ) {
        iterator_type last = std::find( e.first, e.last, '\n' );
        const std::string error_msg =
            BIO_MAKE_STRING(
                "Parsed " << map.size() << " entries.\n"
                << "expected: " << e.what()
                << "got: \"" << std::string( e.first, last ) << '"' << "\n"
                << "on line: " << boost::spirit::get_line( e.first ) << "\n"
                << "column: " << boost::spirit::get_column( begin, e.first ) << "\n"
            );
        throw std::logic_error( error_msg );
    }
}

} // namespace spirit





template< TransData type >
typename DataTraits< type >::entry_t::map_t &
get_deserialise_or_parse(
    const typename DataTraits< type >::entry_t::map_t & map )
{
    namespace fs = boost::filesystem;

    //get rid of const-ness if deserialising or parsing
    typename DataTraits< type >::entry_t::map_t & non_const_map =
        const_cast< typename DataTraits< type >::entry_t::map_t & >( map )
        ;

    if( map.empty() )
    {
        const fs::path serialised_file( DataTraits< type >::get_serialised_binary_file() );

        deserialise_or_init< true >(
            non_const_map,
            serialised_file,
            boost::bind(
                spirit::parse< type >,
                boost::ref( non_const_map ),
                DataTraits< type >::get_biobase_file()
            )
        );
    }

    return non_const_map;
}






BIO_NS_END

#endif //BIO_BIOBASE_PARSE_SPIRIT_H_
