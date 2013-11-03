#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

using namespace boost::spirit;

struct entry {
    int accession_number;
    entry( int i ) : accession_number( i ) { }
};

typedef boost::shared_ptr< entry > entry_ptr;
typedef std::map< int, entry_ptr > map_type;
typedef std::pair< int, entry_ptr > map_entry;


template< typename T >
T &
deref_shared_ptr( boost::shared_ptr< T > & ptr ) {
    return *ptr;
}


/**
 * Parses a biobase table.
 */
template< typename Iterator >
struct map_parser
: qi::grammar< Iterator, map_type() >
{

    // rules
    qi::rule< Iterator, map_type() > start;
    qi::rule< Iterator, map_entry() > map_entry_rule;
    qi::rule< Iterator, entry_ptr() > entry_rule;

    map_parser() : map_parser::base_type( start ) {
        using namespace qi;
        using boost::phoenix::bind;
        using boost::phoenix::construct;
        using boost::phoenix::new_;
        using qi::_1;

        entry_rule = int_[ _val = construct< entry_ptr >( new_< entry >( _1 ) ) ] ;

        map_entry_rule =
                // make each entry into the value type of the map we will insert it into
                entry_rule[
                   bind( &map_entry::first, _val ) = bind( &entry::accession_number, bind( deref_shared_ptr< entry >, _1 ) ),
                   bind( &map_entry::second, _val ) = _1
                ]
            >   lit( "," )
            ;

        start = +map_entry_rule;
    }
};


int main() {
    using namespace boost::spirit;

    typedef std::string::const_iterator iterator;

    const std::string input = "1,2,3,4,";
    iterator iter = input.begin();

    map_parser< std::string::const_iterator > p;

    map_type map;
    const bool result = phrase_parse( iter, input.end(), p, qi::space, map );

    if( result && iter == input.end() ) {
        std::cout << "Parsing successful!\n";
    } else {
        std::cout << "Parsing failed!\n";
    }

    BOOST_FOREACH( map_entry const e, map ) {
        std::cout << e.first << " " << e.second->accession_number << "\n";
    }

    return 0;
}
