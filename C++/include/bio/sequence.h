#ifndef BIO_SEQUENCE_H_
#define BIO_SEQUENCE_H_

#include "bio/defs.h"
#include "bio/hmm_alphabet.h"
#include "bio/random.h"

#include <boost/iterator/filter_iterator.hpp>

#include <string>
#include <list>
#include <vector>
#include <cassert>


BIO_NS_START


typedef std::string seq_t;
typedef std::list<seq_t> SeqList;
extern const seq_t nucleotides;
extern const seq_t iupac_codes;

inline char complement(char c)
{
	switch(c)
	{
	case 'a': case 'A': return 't';
	case 'c': case 'C': return 'g';
	case 'g': case 'G': return 'c';
	case 't': case 'T': return 'a';
	case 'n': case 'N': return 'n';
	default: throw std::logic_error( BIO_MAKE_STRING( "Bad char: " << c ) );
	}
};

template < typename InsIt >
void
reverse_complement(const seq_t & seq, InsIt insert_it)
{
	for (seq_t::const_reverse_iterator i = seq.rbegin();
		seq.rend() != i;
		++i)
	{
		*insert_it++ = complement(*i);
	}
}

struct is_known_nucleotide {
	bool operator()(char c) const {
		return
			'A' == c || 'C' == c || 'G' == c || 'T' == c
			|| 'a' == c || 'c' == c || 'g' == c || 't' == c;
	}
};


enum NucleoCode {
	NUCLEO_A, //A
	NUCLEO_C, //C
	NUCLEO_G, //G
	NUCLEO_T, //T
	NUCLEO_N, //N
	NUM_NUCLEO_CODES
};
typedef std::vector<NucleoCode> NucleoSeq;

inline NucleoCode char_to_nucleo(char c) {
	switch(c) {
		case 'a': case 'A': return NUCLEO_A;
		case 'c': case 'C': return NUCLEO_C;
		case 'g': case 'G': return NUCLEO_G;
		case 't': case 'T': return NUCLEO_T;
		case 'n': case 'N': return NUCLEO_N;
		default: throw std::logic_error( "Bad nucleotide char" );
	}
}
template <>
struct AlphabetTraits<NucleoCode> {
	static size_t get_size() {
		return 4; //we ignore the NUCLEO_N
	}
	static size_t get_index(NucleoCode a) {
		switch (a) {
			case NUCLEO_A: return 0;
			case NUCLEO_C: return 1;
			case NUCLEO_G: return 2;
			case NUCLEO_T: return 3;
			case NUCLEO_N: return 4;
			default: throw std::logic_error( "Unknown alphabet enum" );
		}
	}
	static size_t get_index(char a) {
		switch (a) {
			case 'a': case 'A': return 0;
			case 'c': case 'C': return 1;
			case 'g': case 'G': return 2;
			case 't': case 'T': return 3;
			case 'n': case 'N': return 4;
			default: throw std::logic_error( "Unknown alphabet char" );
		}
	}
	static NucleoCode get_symbol(size_t index) {
		switch(index) {
			case 0: return NUCLEO_A;
			case 1: return NUCLEO_C;
			case 2: return NUCLEO_G;
			case 3: return NUCLEO_T;
			case 4: return NUCLEO_N;
			default: throw std::logic_error( "Unknown alphabet index" );
		}
	}
	static const char * get_string(NucleoCode a) {
		switch(a) {
			case NUCLEO_A: return "A";
			case NUCLEO_C: return "C";
			case NUCLEO_G: return "G";
			case NUCLEO_T: return "T";
			case NUCLEO_N: return "N";
			default: throw std::logic_error( "Unknown alphabet enum" );
		}
	}
	static char get_char(NucleoCode a) {
		switch(a) {
			case NUCLEO_A: return 'A';
			case NUCLEO_C: return 'C';
			case NUCLEO_G: return 'G';
			case NUCLEO_T: return 'T';
			case NUCLEO_N: return 'N';
			default: throw std::logic_error( "Unknown alphabet enum" );
		}
	}
};

inline
std::ostream &
operator<<(std::ostream & os, NucleoCode a)
{
	os << AlphabetTraits<NucleoCode>::get_string(a);
	return os;
}




class IupacCode {
public:
	enum Value {
		IUPAC_A, //A
		IUPAC_C, //C
		IUPAC_G, //G
		IUPAC_T, //T
		IUPAC_N, //A, C, G or T
		IUPAC_R, //A or G
		IUPAC_S, //C or G
		IUPAC_M, //A or C
		IUPAC_W, //A or T
		IUPAC_Y, //C or T
		IUPAC_K, //G or T
		IUPAC_V, //A, C or G
		IUPAC_U, //UNDEFINED! same as 'N'?
		IUPAC_H, //A, C or T
		IUPAC_D, //A, G or T
		IUPAC_B, //C, G or T
		NUM_IUPAC_CODES
	};

	IupacCode(Value v)
		: _value (v)
	{
	}

	IupacCode(char c)
		: _value (char_to_iupac(c))
	{
	}

	static Value char_to_iupac(char c) {
		switch(c)
		{
			case 'a': case 'A': return IUPAC_A;
			case 'c': case 'C': return IUPAC_C;
			case 'g': case 'G': return IUPAC_G;
			case 't': case 'T': return IUPAC_T;
			case 'n': case 'N': case '.': return IUPAC_N;
			case 'r': case 'R': return IUPAC_R;
			case 's': case 'S': return IUPAC_S;
			case 'm': case 'M': return IUPAC_M;
			case 'w': case 'W': return IUPAC_W;
			case 'y': case 'Y': return IUPAC_Y;
			case 'k': case 'K': return IUPAC_K;
			case 'v': case 'V': return IUPAC_V;
			case 'h': case 'H': return IUPAC_H;
			case 'd': case 'D': return IUPAC_D;
			case 'b': case 'B': return IUPAC_B;
			default: throw std::logic_error( BIO_MAKE_STRING( "Bad IUPAC char: " << c ) );
		}
	}

	static char iupac_to_char( Value value ) {
		switch(value)
		{
			case IUPAC_A: return 'a';
			case IUPAC_C: return 'c';
			case IUPAC_G: return 'g';
			case IUPAC_T: return 't';
			case IUPAC_N: return 'n';
			case IUPAC_R: return 'r';
			case IUPAC_S: return 's';
			case IUPAC_M: return 'm';
			case IUPAC_W: return 'w';
			case IUPAC_Y: return 'y';
			case IUPAC_K: return 'k';
			case IUPAC_V: return 'v';
			case IUPAC_H: return 'h';
			case IUPAC_D: return 'd';
			case IUPAC_B: return 'b';
			default: throw std::logic_error( "Bad IUPAC value" );
		}
	}

	operator char() const { return iupac_to_char(_value); }

	Value get_value() const { return _value; }

	float_t get_freq(char c) const {
		return get_freq(char_to_nucleo(c));
	}

	float_t get_freq(NucleoCode code) const {
		return ((float_t) get_num(code)) / get_num_observations();
	}

	int get_num_observations() const {
		return 12;
	}

	int get_num(char c) const {
		return get_num(char_to_nucleo(c));
	}

	int get_num(NucleoCode code) const {
		switch (_value) {
			case IUPAC_A:
				switch(code) {
					case NUCLEO_A: return 12;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_C:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 12;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_G:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 12;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_T:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 12;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_N:
				switch(code) {
					case NUCLEO_A: return 3;
					case NUCLEO_C: return 3;
					case NUCLEO_G: return 3;
					case NUCLEO_T: return 3;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_R:
				switch(code) {
					case NUCLEO_A: return 6;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 6;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_S:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 6;
					case NUCLEO_G: return 6;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_M:
				switch(code) {
					case NUCLEO_A: return 6;
					case NUCLEO_C: return 6;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_W:
				switch(code) {
					case NUCLEO_A: return 6;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 6;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_Y:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 6;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 6;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_K:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 6;
					case NUCLEO_T: return 6;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_V:
				switch(code) {
					case NUCLEO_A: return 4;
					case NUCLEO_C: return 4;
					case NUCLEO_G: return 4;
					case NUCLEO_T: return 0;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_H:
				switch(code) {
					case NUCLEO_A: return 4;
					case NUCLEO_C: return 4;
					case NUCLEO_G: return 0;
					case NUCLEO_T: return 4;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_D:
				switch(code) {
					case NUCLEO_A: return 4;
					case NUCLEO_C: return 0;
					case NUCLEO_G: return 4;
					case NUCLEO_T: return 4;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;

			case IUPAC_B:
				switch(code) {
					case NUCLEO_A: return 0;
					case NUCLEO_C: return 4;
					case NUCLEO_G: return 4;
					case NUCLEO_T: return 4;
					default: throw std::logic_error( "Bad nucleo" );
				}
				break;
			default: throw std::logic_error( "Bad value" );
		}
	}

	int get_max() const {
		switch (_value) {
			case IUPAC_A:
			case IUPAC_C:
			case IUPAC_G:
			case IUPAC_T:
				return 12;

			case IUPAC_N:
				return 3;

			case IUPAC_R:
			case IUPAC_S:
			case IUPAC_M:
			case IUPAC_W:
			case IUPAC_Y:
			case IUPAC_K:
				return 6;

			case IUPAC_V:
			case IUPAC_H:
			case IUPAC_D:
			case IUPAC_B:
				return 4;
			default: throw std::logic_error( "Bad value" );
		}
	}

	int get_min() const {
		switch (_value) {
			case IUPAC_A:
			case IUPAC_C:
			case IUPAC_G:
			case IUPAC_T:
			case IUPAC_R:
			case IUPAC_S:
			case IUPAC_M:
			case IUPAC_W:
			case IUPAC_Y:
			case IUPAC_K:
			case IUPAC_V:
			case IUPAC_H:
			case IUPAC_D:
			case IUPAC_B:
				return 0;

			case IUPAC_N:
				return 3;
			default: throw std::logic_error( "Bad value" );
		}
	}

	static inline IupacCode complement(char c) {
		return complement(IupacCode(c));
	}

	static inline Value complement(IupacCode code) {
		switch (code._value) {
			case IUPAC_A: return IUPAC_T;
			case IUPAC_T: return IUPAC_A;

			case IUPAC_C: return IUPAC_G;
			case IUPAC_G: return IUPAC_C;

			case IUPAC_R: return IUPAC_Y;
			case IUPAC_Y: return IUPAC_R;

			case IUPAC_S: return IUPAC_W;
			case IUPAC_W: return IUPAC_S;

			case IUPAC_M: return IUPAC_K;
			case IUPAC_K: return IUPAC_M;

			case IUPAC_V: return IUPAC_B;
			case IUPAC_B: return IUPAC_V;

			case IUPAC_H: return IUPAC_D;
			case IUPAC_D: return IUPAC_H;

			case IUPAC_N: return IUPAC_N;
			default: throw std::logic_error( "Bad code" );
		}
	}

private:
	Value _value;
};
typedef std::vector<IupacCode> IupacSeq;



template <class InsIt>
void generate_random_nucleotide_seq(InsIt insert_it, size_t length)
{
	for (size_t i = 0; i < length; ++i, ++insert_it) {
		*insert_it = nucleotides[get_uniform_index(4)];
	}
}



BIO_NS_END


#endif //BIO_SEQUENCE_H_
