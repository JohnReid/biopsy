
#ifndef BIO_MARKOV_MODEL_H_
#define BIO_MARKOV_MODEL_H_

#include "bio/defs.h"

#include <boost/multi_array.hpp>
#include <boost/io/ios_state.hpp>
//#include "boost/multi_array/extent_gen.hpp"

#include <map>
#include <iostream>

BIO_NS_START


/**
A simple markov model. Counts occurences of different symbols.
*/
template <unsigned order, typename count_t = unsigned>
struct MarkovModel
{
	//
	// typedefs
	//
	typedef boost::multi_array< count_t, order + 1 > array_t;
	typedef typename array_t::index index_t;
	typedef boost::array< unsigned, order + 1 > indices_t;
	typedef typename boost::array< index_t, order + 1 > shape_t;


	//
	// members
	//
	array_t counts;
	count_t total_count;




	//
	// methods
	//

	/** Construct a markov model. Initialize all counts to 0. */
	MarkovModel(unsigned alphabet_size);

	unsigned get_alphabet_size() const;
	unsigned get_order() const;
	
	/** Add elements in the sequence [begin, end) to the counts. Use alphabet
	to decide on the alphabet index for each value. Returns # elements added. */
	template <typename It, typename SymbolAlphabet>
	unsigned
	add_to_counts(It begin, It end, SymbolAlphabet alphabet)
	{
		BOOST_ASSERT(alphabet.get_alphabet_size() == get_alphabet_size());

		unsigned original_total = total_count;

		//iterators pointing to start and end of current symbol
		typedef SymbolIterator< It, SymbolAlphabet > SymbolIt;
		for (SymbolIt s(begin, end, alphabet);
			end != s.begin;
			s.increment())
		{
			std::vector< unsigned > index;

			//for each character in the current symbol
			for (It c = s.begin;
				s.end != c;
				++c)
			{
				index.push_back( alphabet(*c) );
			}

			++counts(index);
			++total_count;

			BOOST_ASSERT(std::accumulate(counts.data(), counts.data() + counts.num_elements(), count_t(0)) == total_count);
		}

		return unsigned(total_count - original_total);
	}

	/**
	Calls fn for each count.
	*/
	template <typename SymbolAlphabet, typename Fn>
	void
	for_each_count(SymbolAlphabet alphabet, Fn fn)
	{
		indices_t indices;
		std::fill( indices.begin(), indices.end(), 0 );

		for (bool finished = false; ! finished; )
		{
			fn( indices, counts(indices) );

			//increment the indices
			for (unsigned i = order + 1;
				0 != i;
				--i)
			{
				++indices[i - 1];
				if (SymbolAlphabet::get_alphabet_size() == indices[i - 1])
				{
					indices[i - 1] = 0;

					if (1 == i)
					{
						finished = true;
						break;
					}
				}
				else
				{
					break;
				}
			}
		}
	}


	/**
	Calls fn for each count in a sorted order.
	*/
	template <typename SymbolAlphabet, typename Fn>
	void
	for_each_count_sorted(SymbolAlphabet alphabet, Fn fn)
	{
		ForEachSorter sorter;

		for_each_count< SymbolAlphabet, ForEachSorter & >( alphabet, sorter );

		for (typename ForEachSorter::map_t::const_iterator i = sorter.map.begin();
			sorter.map.end() != i;
			++i)
		{
			fn( i->second, i->first );
		}

	}


	template < typename SymbolAlphabet >
	void print( SymbolAlphabet alphabet, bool sorted = true, std::ostream & stream = std::cout )
	{
		stream << "Total : " << total_count << "\n";
		if (sorted)
		{
			for_each_count_sorted( alphabet, Printer< SymbolAlphabet >( alphabet, total_count, stream ) );
		}
		else
		{
			for_each_count( alphabet, Printer< SymbolAlphabet >( alphabet, total_count, stream ) );
		}
	}



	//
	// classes
	//

	/**
	Helps to sort the counts.
	*/
	struct ForEachSorter
	{
		typedef std::multimap< unsigned, indices_t > map_t;
		map_t map;

		void operator()(const indices_t & indices, count_t count)
		{
			map.insert(typename map_t::value_type(count, indices));
		}
	};

	/**
	Prints the counts.
	*/
	template < typename SymbolAlphabet >
	struct Printer
	{
		std::ostream & stream;
		boost::io::ios_all_saver ias;
		SymbolAlphabet alphabet;
		count_t total_count;

		Printer( SymbolAlphabet alphabet, count_t total_count, std::ostream & stream = std::cout )
			: stream(stream)
			, ias(stream)
			, alphabet(alphabet)
			, total_count(total_count)
		{
			stream.fill( ' ' );
			stream.precision( 5 );
		}

		void operator()(const indices_t & indices, count_t count)
		{
			//print the symbol
			for (typename indices_t::const_iterator i = indices.begin();
				indices.end() != i;
				++i)
			{
				std::cout << alphabet.get_char(*i);
			}
			std::cout.setf( std::ios_base::right, std::ios_base::adjustfield );
			std::cout << ", " << std::setw(8) << count;
			std::cout.setf( std::ios_base::left, std::ios_base::adjustfield );
			std::cout << ", " << (0 != total_count ? double(count) / double(total_count) : 0) << "\n";
		}
	};

	/**
	Iterates over the symbols in a sequence.
	*/
	template <typename It, typename SymbolAlphabet>
	struct SymbolIterator
	{
		It begin; /**< The beginning of the current symbol. */
		It end; /**< The end of the current symbol. */
		It seq_end; /**< The end of the whole sequence we are iterating over. */
		SymbolAlphabet alphabet;

		SymbolIterator( It seq_begin, It seq_end, SymbolAlphabet alphabet )
			: begin(seq_begin)
			, seq_end(seq_end)
			, alphabet(alphabet)
		{
			find_next_whole_symbol();
		}

		void find_next_whole_symbol()
		{
			//keep iterating until we find a whole symbol or the end of the sequence
			for ( ;
				seq_end != begin;
				++begin)
			{
				//initialise end pointer
				end = begin;
				unsigned num_chars_in_symbol = 0;

				//iterate until we have required # of chars or found undefined symbol
				for ( ;
					seq_end != end && num_chars_in_symbol != order + 1;
					++num_chars_in_symbol, ++end)
				{
					if (! alphabet.is_defined(*end))
					{
						//start looking again
						num_chars_in_symbol = 0;
						begin = end;
						break;
					}
				}

				//did we find a whole sequence?
				if (num_chars_in_symbol == order + 1)
				{
					//yes - we are done
					break;
				}
			}
		}

		void increment()
		{
			//are we already at the end?
			if (seq_end == end)
			{
				//yes - set the beginning to point at the end and return
				begin = seq_end;
				return;
			}

			//increment the begin iterator
			++begin;

			//is the new end point well defined?
			if (! alphabet.is_defined(*end))
			{
				//no - find the next whole symbol
				find_next_whole_symbol();
			}
			else
			{
				//increment the end iterator to point just past the end character
				++end;
			}
		}
	};

	//
	// static
	//
	static shape_t get_shape(unsigned alphabet_size);
};



template <unsigned order, typename count_t>
MarkovModel<order, count_t>::MarkovModel(unsigned alphabet_size)
: counts(get_shape(alphabet_size))
, total_count(count_t(0))
{
	//initialize counts to 0
	std::fill(counts.data(), counts.data() + counts.num_elements(), count_t(0));
}



template <unsigned order, typename count_t>
unsigned
MarkovModel<order, count_t>::get_order() const
{
	return order;
}



template <unsigned order, typename count_t>
unsigned
MarkovModel<order, count_t>::get_alphabet_size() const
{
	BOOST_ASSERT(0 != counts.num_dimensions());

	return counts.size();
}



template <unsigned order, typename count_t>
typename MarkovModel<order, count_t>::shape_t
MarkovModel<order, count_t>::get_shape(unsigned alphabet_size)
{
	shape_t shape;
	std::fill(shape.begin(), shape.end(), alphabet_size);

	return shape;
}


/**
Converts DNA bases to index [0,3] and back again. Designed for use with MarkovModel class.
*/
struct DnaSymbolAlphabet
{
	static unsigned get_alphabet_size()
	{
		return 4;
	}

	bool is_defined(char c) const
	{
		return operator()(c) < 4;
	}

	unsigned operator()(char c) const
	{
		switch(c)
		{
		case 'a':
		case 'A':
			return 0;
		case 'c':
		case 'C':
			return 1;
		case 'g':
		case 'G':
			return 2;
		case 't':
		case 'T':
			return 3;
		}

		//return out of range for unknown
		return 4;
	}

	char get_char(unsigned idx) const
	{
		switch (idx)
		{
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default:
			throw std::invalid_argument("Out of range index argument in DnaSymbolAlphabet");
		}
	}
};



BIO_NS_END



#endif //BIO_MARKOV_MODEL_H_

