#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <iostream>
#include <string>

using namespace boost::spirit::qi;
using namespace boost::phoenix;
//using namespace std;

struct A { };

std::ostream &
operator<<( std::ostream & os, A const & a ) {
	return os << "an A";
}

typedef std::string::const_iterator iterator;
typedef A arg_type;
typedef rule< iterator, unused_type( arg_type & ) > symbol_rule;
typedef rule< iterator, unused_type() > noarg_rule;

template< typename Val >
noarg_rule
curry( symbol_rule & r, Val & v ) {
	return r( v );
}

int main() {


//	symbol_rule id, one, two;
//	rule<iterator, A() > start;
//	rule<iterator, unused_type( A & ), locals< symbol_rule * > > locals_;
//	symbols< char, symbol_rule * > keyword;

//	char const * input = "";

	int i = 1;
	A a;
	symbol_rule print_arg = eps[ std::cout << val( "Got: " ) << _r1 << "\n" ];
	auto ref_arg = boost::phoenix::ref( a );
	noarg_rule print_actor = curry( print_arg, ref( a ) );
	bind( &curry< decltype( ref_arg ) >, print_arg, ref_arg );
	//noarg_rule lazy_rule = lazy( val( print_actor ) );
	//noarg_rule lazy_rule =
	//noarg_rule lazy_rule = lazy( boost::phoenix::bind( &curry< decltype( ref_arg ) >, print_arg, ref_arg ) );

	std::string str = "one only\none again\ntwo first,second";
    iterator iter = str.begin();
    iterator end = str.end();
    //bool r = parse( iter, end, lazy( val( boost::phoenix::bind( &curry< decltype( ref_arg ) >, print_arg, ref_arg ) ) ) );
    //std::cout << r << "\n";
}

#if 0

/*=============================================================================
    Copyright (c) 2003 Sam Nabialek
    Copyright (c) 2001-2010 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
///////////////////////////////////////////////////////////////////////////////
//
//  The Nabialek trick.
//
//  [ Sam Nabialek; Somewhere, sometime in 2003... ]    spirit1
//  [ JDG November 17, 2009 ]                           spirit2
//  [ JDG January 10, 2010 ]                            Updated to use rule pointers
//                                                      for efficiency.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <iostream>
#include <string>

namespace client
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    struct A { };

    ///////////////////////////////////////////////////////////////////////////////
    //  Our nabialek_trick grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct nabialek_trick : qi::grammar< Iterator, ascii::space_type, A() >
    {
    	typedef qi::rule<Iterator, ascii::space_type, qi::unused_type( A & )> symbol_rule;

        nabialek_trick() : nabialek_trick::base_type(start)
        {
        	using namespace qi;
            using ascii::alnum;
            using qi::_1;

            id = lexeme[*(ascii::alnum | '_')];
            one = id;
            two = id >> ',' >> id;

            keyword.add
                ("one", &one)
                ("two", &two)
                ;

            locals_ = *( keyword[_a = _1] >> lazy( ( *_a )( _r1 ) ) );

            start = locals_( _val );
        }

        symbol_rule id, one, two;
        qi::rule<Iterator, ascii::space_type, A() > start;
        qi::rule<Iterator, ascii::space_type, qi::unused_type( A & ), qi::locals<symbol_rule*> > locals_;
        qi::symbols<char, qi::rule<Iterator, ascii::space_type, qi::unused_type( A & )>*> keyword;
    };
}

///////////////////////////////////////////////////////////////////////////////
//  Main program
///////////////////////////////////////////////////////////////////////////////
int
main()
{
    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;
    typedef client::nabialek_trick<iterator_type> nabialek_trick;

    nabialek_trick g; // Our grammar

    std::string str = "one only\none again\ntwo first,second";
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bool r = phrase_parse(iter, end, g, space);

    if (r && iter == end)
    {
        std::cout << "-------------------------\n";
        std::cout << "Parsing succeeded\n";
        std::cout << "-------------------------\n";
    }
    else
    {
        std::string rest(iter, end);
        std::cout << "-------------------------\n";
        std::cout << "Parsing failed\n";
        std::cout << "stopped at: \": " << rest << "\"\n";
        std::cout << "-------------------------\n";
    }

    return 0;
}


#endif //0
