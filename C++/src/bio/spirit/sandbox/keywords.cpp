#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/adapt_assoc_struct.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <boost/spirit/repository/include/qi_kwd.hpp>
#include <boost/spirit/repository/include/qi_keywords.hpp>
#include <iostream>
#include <string>

struct A {
	int i;
	double f;
};

namespace keys {
struct i_key;
struct f_key;
}

//BOOST_FUSION_ADAPT_STRUCT( A,
//    (int, i)
//    (double, f)
//)

BOOST_FUSION_ADAPT_ASSOC_STRUCT( A,
    (int, i, keys::i_key)
    (double, f, keys::f_key)
)

int main() {
	using namespace boost::spirit::qi;
	using boost::spirit::repository::qi::kwd;
    //using namespace boost::fusion;
    using boost::phoenix::at_key;
    using namespace keys;

    typedef std::string::const_iterator iterator;

	rule< iterator, ascii::blank_type, A() > keywords;
	keywords =
	        kwd( "KEY1" )[ int_[ at_key< i_key >( _val ) = _1 ] ]
        /   kwd( "KEY2" )[ double_[ at_key< f_key >( _val ) = _1 ] ]
        /   kwd( "KEY3" )[ +char_ ]
	    ;

	const std::string input( "KEY1 1 KEY2 2.5 KEY3 string");
	A a;
	iterator it = input.begin();
	if( ! phrase_parse( it, input.end(), keywords, ascii::blank, a ) ) {
        std::cerr << "Could not parse input.\n";
        return -1;
	}
	if( input.end() != it ) {
	    std::cerr << "Did not parse all input.\n";
	    return -2;
	}
	std::cout << "A: " << a.i << ", " << a.f << "\n";

	return 0;
}
