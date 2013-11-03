#include <iostream>
#include <string>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/fusion/sequence.hpp>

using namespace std;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;

using boost::fusion::at_key;

typedef string::const_iterator iter_type;

struct Couple {
    int a;
    int b;
    Couple() : a(0), b(0) {}
};

namespace keys {
    struct first;
    struct second;
}

BOOST_FUSION_ADAPT_ASSOC_STRUCT(
    Couple,
    (int, a, keys::first)
    (int, b, keys::second)
    )


struct G: qi::grammar< iter_type, Couple(), ascii::space_type >
{
    G() : G::base_type( start_rule ) {
        using qi::_val;
        using qi::_1;
        using qi::_2;
        using boost::phoenix::at_key;

        start_rule =
                        ( "first" >> qi::int_
                                [ at_key<keys::first>(_val) = _1 ]
                        )
                    ^
                        ( "second" >> qi::int_
                                [ at_key<keys::second>(_val) = _1 ]
                        );
    }

    qi::rule< iter_type, Couple(), ascii::space_type > start_rule;
};

int main() {
    Couple couple;
    string example = "second 2 first 1";
    iter_type begin( example.begin() );
    iter_type end( example.end() );

    // test at_key -- compiles with no error
    at_key<keys::second>(couple) = 5;

    bool ok = qi::phrase_parse( begin, end, G(), ascii::space, couple );
    if ( ok )
        cout << couple.a << " " << couple.b << endl;
    else
        cout << "Parse failed" << endl;

    return 0;
}
