#define BOOST_SPIRIT_USE_PHOENIX_V3 // Use correct spirit version

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>

using namespace boost::spirit;

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
//        typedef forward_iterator_type iterator_type;
//        iterator_type begin = fwd_begin;
//        iterator_type end = fwd_end;
    typedef pos_iterator_type iterator_type;
    iterator_type begin = pos_begin;
    iterator_type end = pos_end;

		qi::phrase_parse(
			begin,
			end,
			lit(""),
			qi::blank
		);

    return 0;
}
