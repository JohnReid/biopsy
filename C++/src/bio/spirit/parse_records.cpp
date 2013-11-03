#include "bio-pch.h"
#include <bio/spirit/transfac_qi.h>

#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_variable.hpp>

#include <vector>
#include <istream>
#include <sstream>
#include <iostream>

using namespace BIO_NS;
using namespace BIO_NS::spirit;
//using namespace boost::spirit::qi;
//using namespace boost::spirit::ascii;


struct printer
{
    typedef boost::spirit::utf8_string string;

    void element( string const & tag, string const & value, int depth ) const
    {
        for (int i = 0; i < (depth*4); ++i) // indent to depth
            std::cout << ' ';

        std::cout << "tag: " << tag;
        if (value != "")
            std::cout << ", value: " << value;
        std::cout << std::endl;
    }
};

void print_info( boost::spirit::info const & what )
{
    using boost::spirit::basic_info_walker;

    printer pr;
    basic_info_walker<printer> walker(pr, what.tag, 0);
    boost::apply_visitor(walker, what.value);
}

// main function
int main( int, const char ** )
{
	// iterate over stream input
	typedef std::istreambuf_iterator< char > base_iterator_type;
	base_iterator_type in_begin( std::cin );

	// convert input iterator to forward iterator, usable by spirit parser
	typedef boost::spirit::multi_pass< base_iterator_type > forward_iterator_type;
	forward_iterator_type fwd_begin = boost::spirit::make_default_multi_pass( in_begin );
	forward_iterator_type fwd_end;

	// wrap forward iterator with position iterator, to record the position
    typedef boost::spirit::classic::position_iterator2< forward_iterator_type > pos_iterator_type;
    pos_iterator_type pos_begin( fwd_begin, fwd_end, "<stdin>" );
    pos_iterator_type pos_end;

	// define the iterator type we will use.
//		typedef forward_iterator_type iterator_type;
//		iterator_type begin = fwd_begin;
//		iterator_type end = fwd_end;
	typedef pos_iterator_type iterator_type;
	iterator_type begin = pos_begin;
	iterator_type end = pos_end;

	try
	{

		// create some parsers we will use.
		factorlink_parser< iterator_type > factorlink;
		identifier_parser< iterator_type > identifier;
		sequence_parser< iterator_type > sequence;
		dbref_parser< iterator_type > dbref_;

		using namespace qi;
		using qi::_1;
		using qi::ascii::blank;

		qi::phrase_parse(
			begin,
			end,
			+(  (	(lit("BF  ") >> factorlink[print_factorlink])
				|   (lit("ID  ") >> identifier[std::cout << _1 << "\n"])
				|   (lit("SQ  ") >> sequence[std::cout << _1 << "\n"])
				|   (lit("DR  ") >> dbref_[std::cout << _1 << "\n"])
				) >> qi::eol
			),
			blank
		);

		if( end != begin ) {
			throw std::logic_error( "Did not parse all input." );
		}

	} catch ( qi::expectation_failure< iterator_type > const & e ) {
		std::cerr << "expected: "; print_info( e.what_);
		iterator_type last = std::find( e.first, e.last, '\n' );
		std::cerr << "got: \"" << std::string( e.first, last ) << '"' << std::endl;
	} catch( const std::exception & e ) {
		std::cerr << "Exception: " << e.what() << std::endl;
		return -1;
	}
	return 0;
}
