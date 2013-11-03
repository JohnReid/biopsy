
#include <bio/markov_model.h>
#include <bio/sequence.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;

#include <iostream>
using namespace std;

//#define VERBOSE_CHECKING

const SeqList symbol_seqs = boost::assign::list_of
	("")
	("N")
	("A")
	("AN")
	("NA")
	("NAC")
	("AC")
	("ACG")
	("ACGT")
	("NACGT")
	("ANCGT")
	("ACNGT")
	("ACGNT")
	("ACGTN")
	("ACGTNACGT")
	;


	
template <typename MM>
struct MMCounter
{
	MM & mm;

	MMCounter(MM & mm)
		: mm(mm)
	{
	}

	void operator()(const seq_t & sequence)
	{
		mm.add_to_counts(sequence.begin(), sequence.end(), DnaSymbolAlphabet());
	}
};

BOOST_TEST_DONT_PRINT_LOG_VALUE( seq_t::const_iterator );


template <unsigned order>
void
check_symbol_iterator()
{
	cout << "******* check_symbol_iterator() : order = " << order << endl;

	const SeqList symbol_seqs = boost::assign::list_of
		("")
		("N")
		("A")
		("AN")
		("NA")
		("NAC")
		("AC")
		("ACG")
		("ACGT")
		("NACGT")
		("ANCGT")
		("ACNGT")
		("ACGNT")
		("ACGTN")
		("ACGTNACGT")
		;

	for (SeqList::const_iterator s = symbol_seqs.begin();
		symbol_seqs.end() != s;
		++s)
	{
		typedef
			typename MarkovModel<order>::template SymbolIterator<
				typename seq_t::const_iterator,
				DnaSymbolAlphabet>
			SymbolIt;

#ifdef VERBOSE_CHECKING
		cout << "Sequence: \"" << s->c_str() << "\"\n";
#endif

		for (SymbolIt it(s->begin(), s->end(), DnaSymbolAlphabet());
			s->end() != it.begin;
			it.increment())
		{

#ifdef VERBOSE_CHECKING
			copy( it.begin, it.end, ostream_iterator< char >( cout ));
			cout << "\n";
#endif

			BOOST_CHECK_EQUAL(it.end, std::find(it.begin, it.end, 'N'));
		}

#ifdef VERBOSE_CHECKING
		cout << "\n";
#endif

	}
}

template <unsigned order>
void
check_markov_model()
{
	cout << "******* check_markov_model() : order = " << order << endl;

	typedef MarkovModel<order> mm_t;
	mm_t mm(4);
	std::for_each(symbol_seqs.begin(), symbol_seqs.end(), MMCounter<mm_t>(mm));

#ifdef VERBOSE_CHECKING
	typedef boost::array< unsigned, order + 1> indices_t;
	indices_t indices;
	std::fill( indices.begin(), indices.end(), 0 );
	DnaSymbolAlphabet alphabet;

	cout << "Total : " << mm.total_count << "\n";
	while (true)
	{
		//print the symbol
		for (indices_t::const_iterator i = indices.begin();
			indices.end() != i;
			++i)
		{
			cout << alphabet.get_char(*i);
		}
		cout << " : " << mm.counts(indices) << "\n";


		//increment the indices
		for (unsigned i = order + 1;
			0 != i;
			--i)
		{
			++indices[i - 1];
			if (DnaSymbolAlphabet::get_alphabet_size() == indices[i - 1])
			{
				indices[i - 1] = 0;

				if (1 == i)
				{
					return;
				}
			}
			else
			{
				break;
			}
		}
	}
#endif
}

void
register_markov_model_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_symbol_iterator<0> ), 0);
	test->add( BOOST_TEST_CASE( &check_symbol_iterator<1> ), 0);
	test->add( BOOST_TEST_CASE( &check_symbol_iterator<2> ), 0);
	test->add( BOOST_TEST_CASE( &check_markov_model<0> ), 0);
	test->add( BOOST_TEST_CASE( &check_markov_model<1> ), 0);
	test->add( BOOST_TEST_CASE( &check_markov_model<2> ), 0);
}


