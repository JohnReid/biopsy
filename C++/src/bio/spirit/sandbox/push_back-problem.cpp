#define BOOST_SPIRIT_USE_PHOENIX_V3

#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/define_struct.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <vector>

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

BOOST_FUSION_DEFINE_STRUCT(
    ,
    my_struct,
    (std::vector< int >, member)
)

int
main( int argc, char * argv[] )
{
    qi::rule< const char *, my_struct() > start = qi::int_[ phoenix::push_back( phoenix::bind( &my_struct::member, qi::_val ), qi::_1 ) ];

    return 0;
}
