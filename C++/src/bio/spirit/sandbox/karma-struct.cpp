#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>
#include <iostream>

//#define HAS_PTR_CONSTRUCTOR

struct number_container {
	number_container( float a, float b, float c )
	: a(a), b(b), c(c) { }

	number_container()
	: a(0.), b(0.), c(0.) { }

#ifdef HAS_PTR_CONSTRUCTOR
	number_container( float * p ) {
		a = p[0];
		b = p[1];
		c = p[2];
	}
#endif //HAS_PTR_CONSTRUCTOR

	float a, b, c;
};

BOOST_FUSION_ADAPT_ADT(
    number_container,
    (bool, bool, true, /**/)
    (int, int, 1, /**/)
//    (char, char, 'x', /**/)
//    (float, float, 0, /**/)
//    (double, double, 1, /**/)
//    (std::string, std::string, "string", /**/)
//    (float, float, obj.a, /**/)
//    (float, float, obj.b, /**/)
//    (float, float, obj.c, /**/)
)

namespace karma = boost::spirit::karma;

template< typename Iterator >
struct test_generator
: karma::grammar< Iterator, number_container() >
{
    karma::rule< Iterator, number_container() > start;

    test_generator() : test_generator::base_type( start )
    {
        using namespace karma;

        start =
        	eps
        	<< bool_
        	<< int_
//        	char_
//        	<< float_
//        	<< double_
//        	<< *char_
//        	float_ << lit("  ")
//        	<< float_ << lit("  ")
//        	<< float_ << lit("  ")
        	;
    }
};

int main() {

#ifdef HAS_PTR_CONSTRUCTOR
	std::cout << "Has ptr constructor.\n";
	float p[3] = {0., 1., 2.};
//	p[0] = 0.;
//	p[1] = 1.;
//	p[2] = 2.;
	number_container x(p);
#else
	std::cout << "Does not have ptr constructor.\n";
	number_container x(0., 1., 2.);
#endif //HAS_PTR_CONSTRUCTOR

	typedef std::ostream_iterator< char > output_iterator;
	test_generator< output_iterator > generator;
	output_iterator sink( std::cout );
	if( ! karma::generate(
			sink,
			generator,
			x
		)
	) {
		std::cerr << "Could not generate output.\n";
	}

	return 0;
}
