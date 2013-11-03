// #define BOOST_SPIRIT_DEBUG
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <map>

namespace qi    = boost::spirit::qi;
namespace karma = boost::spirit::karma;
namespace phx   = boost::phoenix;

typedef std::map<std::string, std::string> pairs_t;

template <typename It, typename Skipper = qi::space_type>
    struct parser : qi::grammar<It, pairs_t(), Skipper>
{
    parser() : parser::base_type(start)
    {
        using namespace qi;

        name  = lexeme [ +(graph - '_') ]; // or just char_("ABCDEFG") if it's that simple
        entry = lexeme [ raw [ +(name >> '_') ] >> name ];
        start = *entry;
        BOOST_SPIRIT_DEBUG_NODE(name);
        BOOST_SPIRIT_DEBUG_NODE(entry);
        BOOST_SPIRIT_DEBUG_NODE(start);
    }

  private:
    qi::rule<It, std::string()> name;
    qi::rule<It, std::pair<std::string, std::string>(), Skipper> entry;
    qi::rule<It, pairs_t(), Skipper> start;
};

template <typename C, typename Skipper>
    bool doParse(const C& input, const Skipper& skipper)
{
    auto f(std::begin(input)), l(std::end(input));

    parser<decltype(f), Skipper> p;
    pairs_t data;

    try
    {
        bool ok = qi::phrase_parse(f,l,p,skipper,data);
        if (ok)
        {
            std::cout << "parse success\n";
            std::cout << "data: " << karma::format_delimited(karma::auto_ % ';', ' ', data) << '\n';
        }
        else    std::cerr << "parse failed: '" << std::string(f,l) << "'\n";

        if (f!=l) std::cerr << "trailing unparsed: '" << std::string(f,l) << "'\n";
        return ok;
    } catch(const qi::expectation_failure<decltype(f)>& e)
    {
        std::string frag(e.first, e.last);
        std::cerr << e.what() << "'" << frag << "'\n";
    }

    return false;
}

int main()
{
    const std::string input = "A_B A_B_C A_B_C_D A_B_C_D_E";
    bool ok = doParse(input, qi::blank);

    return ok? 0 : 255;
}

