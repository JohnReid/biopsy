/**
 * Minimal example to show invalid initialisation that happens with boost 1.52 and above.
 * I think the compilation error is due to boost changing their implementation of result_of.
 */

#define BOOST_SPIRIT_USE_PHOENIX_V3

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

using namespace boost::spirit;

struct entry {
    int accession_number;
};

typedef boost::shared_ptr< entry > entry_ptr;
typedef std::map< int, entry_ptr > map_type;
typedef std::pair< int, entry_ptr > map_entry;


template< typename T >
T &
deref_shared_ptr( boost::shared_ptr< T > & ptr ) {
    return *ptr;
}


template< typename Iterator >
struct map_grammar
: qi::grammar< Iterator, map_type() >
{
    // rules
    qi::rule< Iterator, map_type() > start;
    qi::rule< Iterator, map_entry() > map_entry_rule;
    qi::rule< Iterator, entry_ptr() > entry_rule;

    map_grammar() : map_grammar::base_type( start ) {
        using namespace qi;
        using boost::phoenix::bind;
        using qi::_1;

        map_entry_rule =
                // make each entry into the value type of the map we will insert it into
                entry_rule[
                   bind( &map_entry::first, _val ) = bind( &entry::accession_number, *_1 )
                ]
            ;
    }
};


int main() {
    map_grammar< const char * > p;
}
