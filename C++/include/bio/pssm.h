

#ifndef BIO_PSSM_H_
#define BIO_PSSM_H_

#include "bio/defs.h"
#include "bio/sequence.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#include <vector>


BIO_NS_START




class PssmEntry
{
public:
	/** counts is an array of length 4 specifying the counts of the nucleotides. */
	PssmEntry(float_t * counts);

	PssmEntry(
		float_t num_a,
		float_t num_c,
		float_t num_g,
		float_t num_t);

	/** Initialise with zero counts. */
	PssmEntry();

	float_t get_num_observations() const;
	float_t get_max() const;
	float_t get_min() const;
	float_t get_freq( char c, float_t pseudo_count = 0.0 ) const;
	float_t get_conservation_information() const;
	inline float_t get_count(char c) const { return counts[index_of(c)]; }
	inline float_t get_score(char c) const { return scores[index_of(c)]; }

protected:
	float_t counts[4]; //counts of how many nucleotides in this entry
	float_t scores[4]; //score for each nucleotide

	static inline size_t index_of(char c)
	{
		switch (c)
		{
			case 'a': case 'A': return 0;
			case 'c': case 'C': return 1;
			case 'g': case 'G': return 2;
			case 't': case 'T': return 3;
		}
		throw std::logic_error( BIO_MAKE_STRING( "Bad char: " << c ) );
	}

	void init();

private:
    friend class boost::serialization::access;
	/** De/serialise from/to archive. */
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		for (int i = 0; 4 != i; ++i)
		{
	        ar & counts[i];
	        ar & scores[i];
		}
    }		 
};

class Pssm : public std::vector<PssmEntry>
{
public:

	typedef boost::shared_ptr<Pssm> ptr_t;
	typedef float_t score_t;

	Pssm();
	~Pssm();

	score_t
	score(
		seq_t::const_iterator seq_begin,
		bool match_complement = false) const;

private:
    friend class boost::serialization::access;
	/** De/serialise from/to archive. */
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
        ar & boost::serialization::base_object<std::vector <PssmEntry> >(*this);
    }		 
};



struct PssmScorer
	: std::binary_function< seq_t::const_iterator, bool, float_t >
{
	const Pssm * pssm;

	PssmScorer( const Pssm & pssm )
		: pssm( boost::addressof( pssm ) )
	{
	}

	PssmScorer()
		: pssm( 0 )
	{
	}

	float_t operator()(
		seq_t::const_iterator seq_begin,
		bool match_complement) const
	{
		check_pointer();
		return pssm->score( seq_begin, match_complement );
	}

	unsigned get_num_bases() const
	{
		check_pointer();
		return pssm->size();
	}

private:
	void check_pointer() const
	{
		if( 0 == pssm )
		{
			throw std::logic_error( "Null pointer in PssmScorer" );
		}
	}
};


inline
std::ostream &
operator<<(std::ostream & os, PssmEntry const& pssm_entry)
{
	os
		<< pssm_entry.get_count('a') << "\t"
		<< pssm_entry.get_count('c') << "\t"
		<< pssm_entry.get_count('g') << "\t"
		<< pssm_entry.get_count('t') << "\t"
		<< "\t"
		<< pssm_entry.get_score('a') << "\t"
		<< pssm_entry.get_score('c') << "\t"
		<< pssm_entry.get_score('g') << "\t"
		<< pssm_entry.get_score('t');
	return os;
}

inline
std::ostream &
operator<<(std::ostream & os, Pssm const& pssm)
{
	std::copy(pssm.begin(), pssm.end(), std::ostream_iterator< PssmEntry >(os, "\n"));
	return os;
}


BIO_NS_END


#endif //BIO_PSSM_H_


