#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <iostream>
#include <vector>
#include <string>

namespace qi = boost::spirit::qi;


template< typename P >
void
test_parser(
    char const * input,
    P const & p,
    bool full_match = true
) {
	using namespace ::boost::spirit::qi;
    char const * f( input );
    char const * l( f + strlen( f ) );
    if( !  parse( f, l, p ) ) {
    	std::cerr << "Problem parsing \"" << input << "\"" << "\n";
    } else {
    	std::cerr << "Parsed \"" << input << "\"" << "\n";
    }
    if( full_match && f != l ) {
    	std::cerr << "Did not fully match \n    \"" << input << "\""
    		<< "\n  only matched \n    \"" << std::string( input, f ) << "\"" << "\n";
    }
}

template< typename P, typename T >
T &
test_phrase_parser_attr(
    char const * input,
    P const & p,
    T & attr,
    bool full_match = true
) {
	using namespace ::boost::spirit::qi;
    char const * f( input );
    char const * l( f + strlen( f ) );
    if( ! phrase_parse( f, l, p, ascii::space, attr ) ) {
    	std::cerr << "Problem parsing \"" << input << "\"" << "\n";
    }
    if( full_match && f != l ) {
    	std::cerr << "Did not fully match \"" << input << "\""
    		<< " only matched \"" << std::string( input, f ) << "\"" << "\n";
    }
    return attr;
}

struct simple {
	std::string a;
	std::string b;
	std::string c;
};

BOOST_FUSION_ADAPT_STRUCT(
    simple,
    (std::string, a)
    (std::string, b)
    (std::string, c)
)

template< typename Iterator >
struct sandbox_parser
: qi::grammar< Iterator, simple() >
{
	sandbox_parser() : sandbox_parser::base_type( start )
    {
        using namespace qi;

        start =
        	(string("A") >> string("B") >> string("C"))
        	|
        	(attr("A") >> string("$") >> attr("C"))
        	;
    }

    qi::rule< Iterator, simple() > start;
};

template < typename Iterator >
struct recursive_parser : qi::grammar< Iterator >
{
	qi::rule< Iterator > start;
	qi::rule< Iterator > end;
	recursive_parser() : recursive_parser::base_type( start )
	{
		using namespace qi;
		end = string("_end") | start;
		start = +(char_ - '_') >> end;
	}
};

void
test_recursive() {
	test_parser( "a_end", recursive_parser< char const * >() );
	test_parser( "a_b_end", recursive_parser< char const * >() );
}

void
test_simple() {
	{
		simple attr;
		sandbox_parser< char const * > parser;
		test_phrase_parser_attr(
			"ABC",
			parser,
			attr
		);
		std::cout << attr.a << attr.b << attr.c << "\n";
	}

	{
		simple attr;
		sandbox_parser< char const * > parser;
		test_phrase_parser_attr(
			"$",
			parser,
			attr
		);
		std::cout << attr.a << attr.b << attr.c << "\n";
	}
}

int
main( int argc, char const *argv[] ) {
	test_recursive();

	return 0;
}
