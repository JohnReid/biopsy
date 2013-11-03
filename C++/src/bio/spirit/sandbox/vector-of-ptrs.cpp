#include <boost/spirit/include/karma.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

typedef boost::shared_ptr< int > ptr;
typedef std::vector< ptr > vec;

namespace boost {
namespace spirit {
namespace traits {

// specialise how iterators into containers of pointers are dereferenced
template<>
struct deref_iterator< karma::detail::indirect_iterator< vec::const_iterator > >
{
	typedef karma::detail::indirect_iterator< vec::const_iterator > It;
	typedef int type;

	static type call( It const & it ) {
		return **it;
	}
};

} // namespace traits
} // namespace spirit
} // namespace boost



int
main() {
	vec v;
	v.push_back( ptr( new int( 1 ) ) );
	v.push_back( ptr( new int( 2 ) ) );
	v.push_back( ptr( new int( 3 ) ) );

	using namespace boost::spirit::karma;
	using boost::spirit::ascii::space;

	typedef std::ostream_iterator< char > iterator;
	rule< iterator, int() > int_line = lit("INT=") << int_ << ',';
	generate_delimited(
		iterator( std::cout ),
		*int_line << eol,
		// *( lit("INT: ") << int_ << eol ) << eol, // this does not compile
		space,
		v
	);

	return 0;
}
