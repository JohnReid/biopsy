#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>
#include <iostream>

struct empty_struct {
	empty_struct() { }
};

BOOST_FUSION_ADAPT_ADT(
    empty_struct,
    (bool, bool, true, /**/)
    //(bool, bool, true, /**/) // this compiles
    (int, int, 1, /**/) // this doesn't
)

namespace karma = boost::spirit::karma;

template< typename Iterator >
struct empty_generator
: karma::grammar< Iterator, empty_struct() >
{
    karma::rule< Iterator, empty_struct() > start;

    empty_generator() : empty_generator::base_type( start )
    {
        using namespace karma;
        //start = bool_ << bool_; // this compiles
        start = bool_ << int_; // this doesn't
    }
};

int main() {

	typedef std::ostream_iterator< char > output_iterator;
	empty_generator< output_iterator > generator;
	output_iterator sink( std::cout );
	karma::generate(
		sink,
		generator,
		empty_struct()
	);

	return 0;
}
