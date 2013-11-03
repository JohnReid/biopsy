#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
#include <boost/phoenix/function/adapt_function.hpp>

namespace qi      = boost::spirit::qi;
namespace phx     = boost::phoenix;
namespace classic = boost::spirit::classic;
namespace ascii   = boost::spirit::ascii;

typedef std::vector<double> data_t;

///////// USING A FREE FUNCTION
//
template <typename Grammar, typename Range>
    double doStuff_(Grammar &grammar, Range pos_range)
{
    // for efficiency, cache adhoc grammar:
    static const qi::rule   <typename Range::iterator, double()> r_double = qi::double_;
    static const qi::grammar<typename Range::iterator, double()> g_double(r_double); // caching just the rule may be enough, actually

    double value = 0;
    qi::parse(pos_range.begin(), pos_range.end(), g_double, value);

    std::cout << "debug ('" << grammar.name() << "') at "
       << pos_range.begin().get_position().file   << ":"
       << pos_range.begin().get_position().line   << ":["
       << pos_range.begin().get_position().column << ".."
       << pos_range.end  ().get_position().column << "]\t"
       << "'" << std::string(pos_range.begin(),pos_range.end()) << "'\t = "
       << value
       << '\n';

    return value;
}

BOOST_PHOENIX_ADAPT_FUNCTION(double, doStuff, doStuff_, 2)

template <typename Iterator, typename Skipper>
struct parse_grammar : qi::grammar<Iterator, data_t(), Skipper>
{
    parse_grammar()
        : parse_grammar::base_type(start_p, "start_p")
    {
        using qi::raw;
        using qi::double_;
        using qi::_1;
        using qi::_val;
        using qi::eoi;
        using phx::push_back;

        value_p = raw [ double_ ] [ _val = doStuff(phx::ref(*this), _1) ];
        start_p = value_p % ',' > eoi;

        // // To use without the semantic action (more efficient):
        // start_p = double_ % ',' >> eoi;
    }

    qi::rule<Iterator, data_t::value_type(), Skipper> value_p;
    qi::rule<Iterator, data_t(), Skipper> start_p;
};

// implementation
data_t parse(std::istream& input, const std::string& filename)
{
    // iterate over stream input
    typedef std::istreambuf_iterator<char> base_iterator_type;
    base_iterator_type in_begin(input);

    // convert input iterator to forward iterator, usable by spirit parser
    typedef boost::spirit::multi_pass<base_iterator_type> forward_iterator_type;
    forward_iterator_type fwd_begin = boost::spirit::make_default_multi_pass(in_begin);
    forward_iterator_type fwd_end;

    // wrap forward iterator with position iterator, to record the position
    typedef classic::position_iterator2<forward_iterator_type> pos_iterator_type;
    pos_iterator_type position_begin(fwd_begin, fwd_end, filename);
    pos_iterator_type position_end;

    parse_grammar<pos_iterator_type, ascii::space_type> gram;

    data_t output;
    // parse
    try
    {
        if (!qi::phrase_parse(
                position_begin, position_end,  // iterators over input
                gram,                          // recognize list of doubles
                ascii::space,                  // comment skipper
                output)                        // <-- attribute reference
           )
        {
            std::cerr << "Parse failed at "
               << position_begin.get_position().file   << ":"
               << position_begin.get_position().line   << ":"
               << position_begin.get_position().column << "\n";
        }
    }
    catch(const qi::expectation_failure<pos_iterator_type>& e)
    {
        const classic::file_position_base<std::string>& pos = e.first.get_position();
        std::stringstream msg;
        msg << "parse error at file " << pos.file
            << " line "               << pos.line
            << " column "             << pos.column
            << "\n\t'"                << e.first.get_currentline()
            << "'\n\t "               << std::string(pos.column, ' ') << "^-- here";

        throw std::runtime_error(msg.str());
    }

    return output;
}

int main()
{
    std::istringstream iss(
            "1, -3.4 ,3.1415926\n"
            ",+inF,-NaN  ,\n"
            "2,-.4,4.14e7\n");

    data_t parsed = parse(iss, "<inline-test>");

    std::cout << "Done, parsed " << parsed.size() << " values ("
        << "min: " << *std::min_element(parsed.begin(), parsed.end()) << ", "
        << "max: " << *std::max_element(parsed.begin(), parsed.end()) << ")\n";
}
