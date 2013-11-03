#include <boost/spirit/include/qi.hpp>

int main() {
	using namespace boost::spirit;

	qi::rule< char const *, ascii::blank_type > blank_skipper;
	qi::rule< char const * > implicit_lexeme;
	qi::rule< char const *, ascii::blank_type > blank_compound = blank_skipper | implicit_lexeme;
	qi::rule< char const * > lexeme_compound = qi::skip( ascii::blank )[ blank_skipper ] | implicit_lexeme;

	return 0;
}
