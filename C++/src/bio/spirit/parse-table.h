#include <boost/config/warning_disable.hpp>

#include <bio/spirit/transfac_qi.h>
#include <bio/spirit/transfac_karma.h>

#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>
#include <boost/spirit/home/classic/iterator/position_iterator.hpp>

#include <vector>
#include <istream>
#include <sstream>
#include <iostream>


using namespace BIO_NS;
using namespace BIO_NS::spirit;

struct printer
{
    typedef boost::spirit::utf8_string string;

    void element( string const & tag, string const & value, int depth ) const
    {
        for (int i = 0; i < (depth*4); ++i) // indent to depth
            std::cout << ' ';

        std::cout << "tag: " << tag;
        if (value != "")
            std::cout << ", value: " << value;
        std::cout << std::endl;
    }
};

void print_info( boost::spirit::info const & what )
{
    using boost::spirit::basic_info_walker;

    printer pr;
    basic_info_walker<printer> walker(pr, what.tag, 0);
    boost::apply_visitor(walker, what.value);
}

// main function
template< typename Entry >
int
parse_table( int argc, char const * argv[] )
{
    if( argc < 2 ) {
        std::cerr << "USAGE: " << argv[ 0 ] << " <TABLE FILE>\n";
        return -4;
    }
    std::ifstream input( argv[1] );
    if( ! input ) {
        std::cerr << "Could not open \"" << argv[ 1 ] << "\"\n";
        return -3;
    }
    input.unsetf( std::ios::skipws ); // disable skipping of whitespace

    // iterate over stream input using boost::spirit forward iterator for streams
    typedef boost::spirit::istream_iterator base_iterator_type;
    base_iterator_type input_begin( input );

    // wrap stream iterator with position iterator, to record the position
    typedef boost::spirit::line_pos_iterator< base_iterator_type > iterator_type;
    iterator_type begin( input_begin );
    iterator_type end;

    typename Entry::map_t entries;
    try
    {
        typedef typename get_entry_parser< Entry, iterator_type >::type parser;
        biobasetable_parser< iterator_type, parser > table;

        using namespace qi;
        using qi::ascii::blank;

        iterator_type i = begin;
        qi::parse(
            i,
            end,
            table,
            entries
        );
        std::cout << "Parsed " << entries.size() << " entries.\n";
        typedef std::ostream_iterator< char > output_iterator;
        typedef typename get_entry_generator< Entry, output_iterator >::type generator_type;
        generator_type generator;
        BOOST_FOREACH( typename Entry::map_t::value_type const & entry, entries ) {
            if( ! karma::generate(
                    output_iterator( std::cout ),          // destination: output iterator
                    generator,                             // the generator
                    *(entry.second)                        // the data to output
                )
            ) {
                throw std::logic_error( "Could not generate output for entry." );
            }
//            BOOST_FOREACH( PssmEntry const & entry, matrix.pssm ) {
//                std::cout
//                    << entry.get_count( 'a' ) << "  "
//                    << entry.get_count( 'c' ) << "  "
//                    << entry.get_count( 'G' ) << "  "
//                    << entry.get_count( 'T' ) << "  "
//                    << "\n";
//            }
        }

        if( end != i ) {
            throw std::logic_error(
                BIO_MAKE_STRING( "Did not parse all input: Got to line "
                    << boost::spirit::get_line( i )
                    << ", column " << boost::spirit::get_column( begin, i ) ) );
        }

    } catch ( qi::expectation_failure< iterator_type > const & e ) {
        iterator_type last = std::find( e.first, e.last, '\n' );
        std::cout << "Parsed " << entries.size() << " entries.\n";
        std::cerr << "expected: "; print_info( e.what_);
        std::cerr << "got: \"" << std::string( e.first, last ) << '"' << "\n";
        std::cerr << "on line: " << boost::spirit::get_line( e.first ) << "\n";
        std::cerr << "column: " << boost::spirit::get_column( begin, e.first ) << "\n";
        return -2;
    } catch( const std::exception & e ) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return -1;
    }
    return 0;
}
