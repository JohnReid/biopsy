#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/adapt_struct_named.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>
#include <iostream>

struct struct_of_floats {
	struct_of_floats() : a(0), b(1) { }
	float a, b;
};

BOOST_FUSION_ADAPT_ADT(
    struct_of_floats,
    (float, float, 0., /**/)
    (float, float, 1., /**/)
)

BOOST_FUSION_ADAPT_STRUCT_NAMED(
    struct_of_floats const, adapted_struct_of_floats,
    (float, a)
    (float, b)
)

int main() {

	using namespace boost::spirit::karma;
	typedef std::ostream_iterator< char > output_iterator;
	output_iterator sink( std::cout );
	rule< output_iterator, struct_of_floats() > rule = float_ << float_;
	// NEXT LINE WILL NOT COMPILE
	//rule< output_iterator, adapted_struct_of_floats() > adapted_rule = float_ << float_;
	generate(
		sink,
		rule,
		struct_of_floats()
	);

	return 0;
}
