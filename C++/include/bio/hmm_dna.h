#ifndef BIO_HMM_DNA_H_
#define BIO_HMM_DNA_H_

#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/background_models.h"
#include "bio/hidden_markov_model.h"
#include "bio/hmm_gen_sequence.h"
#include "bio/hmm_baum_welch.h"
#include "bio/singleton.h"

#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <boost/serialization/access.hpp>

#include <string>
#include <map>



BIO_NS_START

//forward decl
struct ConversionChecker;


/** The alphabet for HMM emissions. */
template< unsigned order >
struct DnaHmmAlphabet
{
	unsigned value;
	DnaHmmAlphabet(unsigned value) : value(value) { }
};




/**
A model of DNA implemented by a HMM of the given order. Will not work with order > 15 as
the unsigned type is not large enough to hold that number of bases.
*/
template< unsigned order >
class DnaHmm
	: public DnaModel
{
	friend struct ConversionChecker;


	/**
	Types
	*/

public:
	/** Functor for use in baum welch. */
	template <bool value>
	struct ConstPredicate
	{
		template <class Type>
		bool operator()(const Type & type) const { return value; }
	};
	typedef DnaHmmAlphabet< order > alphabet_t;
	typedef MarkovState<alphabet_t> state_t;
	typedef HiddenMarkovModel<alphabet_t> model_t;
	typedef HmmSequenceGenerator<alphabet_t> seq_gen_t;
	typedef std::vector<alphabet_t> emission_seq_t;
	typedef std::list<emission_seq_t> emission_seq_list_t;




	/**
	Data
	*/

protected:
	model_t model;






	/**
	Methods
	*/

public:
	/** Construct with given number of states. */
	DnaHmm(unsigned num_states = 1);

	virtual ~DnaHmm();

	/** Train on one sequence. */
	virtual void train(const seq_t & sequence);

	/** Return the per base likelihood under this model of the given sequence. */
	virtual prob_t get_likelihood(const seq_t & sequence);

	/** Train on multiple sequences concurrently. */
	virtual void train(const SequenceCollection & sequences);

	/** Return the per base likelihood under this model of the given sequences. */
	virtual prob_t get_likelihood(const SequenceCollection & sequences);

	/** Generate a random sequence from this model and append it to the sequence. */
	virtual void append_random_sequence(seq_t & seq, unsigned seq_length) const;

    friend class boost::serialization::access;
	/** Serialize. */
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
	{
        // explicitly register base/derived relationship as base class has no serialize method to call
        boost::serialization::void_cast_register<DnaHmm< order >, DnaModel>(0, 0);

		ar & model;
	}

protected:
	static void convert_to_emission(const seq_t & seq, emission_seq_t & emission_seq);
	static void convert_to_emission(const SequenceCollection & seq_list, emission_seq_list_t & emission_seq_list);
	static void convert_to_dna(const emission_seq_t & emission_seq, seq_t & seq);
};


/** Returns x to the power of y. */
inline
unsigned unsigned_power(unsigned x, unsigned y)
{
	unsigned result = 1;
	while (y != 0)
	{
		result *= x;
		--y;
	}
	return result;
}

/** Specialise the alphabet on the order. */
template< unsigned order >
struct AlphabetTraits< DnaHmmAlphabet< order > >
{
	static size_t get_size() { return unsigned_power(4, order + 1); }
	static size_t get_index(DnaHmmAlphabet< order > a) { return a.value; }
	static DnaHmmAlphabet< order > get_symbol(size_t index) { return index; }
};

template< unsigned order >
DnaHmm< order >::DnaHmm(unsigned num_states)
: model(num_states)
{
	if (order > 15)
	{
		throw std::logic_error( "Orders > 15 not supported due to size of unsigned type." );
	}
}

template< unsigned order >
DnaHmm< order >::~DnaHmm()
{ }

template< unsigned order >
void
DnaHmm< order >::train(const seq_t & sequence)
{
	emission_seq_t emission_seq;
	convert_to_emission(sequence, emission_seq);

	baum_welch_single(
		model,
		emission_seq.begin(),
		emission_seq.end(),
		emission_seq.rbegin(),
		emission_seq.rend());

	BOOST_ASSERT(model.is_consistent());
}

template< unsigned order >
prob_t
DnaHmm< order >::get_likelihood(const seq_t & sequence)
{
	emission_seq_t emission_seq;
	convert_to_emission(sequence, emission_seq);

	const prob_t log_prob =
		ForwardBackwardAlgorithm<true>()
			.forward(
				model,
				emission_seq.begin(),
				emission_seq.end()).get_log_probability();

	return
		(0 == emission_seq.size())
			? 1.0
			: std::exp(log_prob / (emission_seq.size() * (order + 1)));
}

template< unsigned order >
void
DnaHmm< order >::train(const SequenceCollection & sequences)
{
	emission_seq_list_t emission_seq_list;
	convert_to_emission(sequences, emission_seq_list);

	BaumWelchMultipleAlgorithm<true> baum_welch;
	baum_welch.run(
		model,
		emission_seq_list.begin(),
		emission_seq_list.end());
}

template< unsigned order >
prob_t
DnaHmm< order >::get_likelihood(const SequenceCollection & sequences)
{
	prob_t prob_sum = 0;
	for (unsigned i = 0; sequences.num_sequences() != i; ++i)
	{
		prob_sum += get_likelihood(sequences.get_sequence(i));
	}
	return (0 == sequences.num_sequences()) ? 1.0 : prob_sum / sequences.num_sequences();
}

template< unsigned order >
void
DnaHmm< order >::append_random_sequence(seq_t & seq, unsigned seq_length) const
{
	seq_gen_t seq_gen(&model);
	emission_seq_t emission_seq;
	while (emission_seq.size() * (order + 1) < seq_length)
	{
		emission_seq.push_back(seq_gen());
	}

	//prepare the sequence
	seq.clear();
	seq.reserve(seq_length);
	convert_to_dna(emission_seq, seq);
	if (seq.size() > seq_length)
	{
		seq.resize(seq_length);
	}

	BOOST_ASSERT(seq.size() == seq_length);
}

template< unsigned order >
void
DnaHmm< order >::convert_to_emission(const seq_t & seq, emission_seq_t & emission_seq)
{
	emission_seq.clear();
	seq_t::const_iterator begin = seq.begin();
	while (begin != seq.end())
	{
		if (! is_known_nucleotide()(*begin))
		{
			throw std::logic_error( "Cannot model unknown bases" );
		}

		alphabet_t next_emission = 0;
		unsigned power = 1;
		unsigned i = 0;
		seq_t::const_iterator end = begin;
		for ( ; end != seq.end() && (order + 1) != i; ++i, ++end)
		{
			next_emission.value += power * AlphabetTraits<NucleoCode>::get_index(*end);
			power *= 4;
		}

		//did we get enough bases to output something?
		if ((order + 1) == i)
		{
			emission_seq.push_back(next_emission);
		}
		else
		{
			BOOST_ASSERT(end == seq.end());
		}

		if (end == seq.end())
		{
			break;
		}

		begin = end;
	}

	BOOST_ASSERT(emission_seq.size() == seq.size() / (order + 1));
}

template< unsigned order >
void
DnaHmm< order >::convert_to_dna(const emission_seq_t & emission_seq, seq_t & seq)
{
	seq.clear();
	for (typename emission_seq_t::const_iterator e = emission_seq.begin();
		emission_seq.end() != e;
		++e)
	{
		alphabet_t emission = *e;
		for (unsigned i = 0; order + 1 != i; ++i)
		{
			seq.push_back(
				AlphabetTraits<NucleoCode>::get_char(
					AlphabetTraits<NucleoCode>::get_symbol(emission.value % 4)));
			emission.value /= 4;
		}
	}
}

template< unsigned order >
void
DnaHmm< order >::convert_to_emission(const SequenceCollection & sequences, emission_seq_list_t & emission_seq_list)
{
	emission_seq_list.clear();
	for (unsigned i = 0; sequences.num_sequences() != i; ++i)
	{
		emission_seq_list.push_back(emission_seq_t());
		convert_to_emission(sequences.get_sequence(i), *(emission_seq_list.rbegin()));
	}
}

DnaModel::ptr_t
create_dna_model(unsigned num_states, unsigned order);




















/**
A collection of hidden markov models of different orders and number of states
*/
struct DnaHmmOrderNumStateMap
	: Singleton< DnaHmmOrderNumStateMap >
{


	/*
	Types
	*/
	typedef DnaModel model_t;
	typedef DnaModel::ptr_t model_ptr_t;
	struct Index
	{
		unsigned num_states;
		unsigned order;
		Index(unsigned num_states = 1, unsigned order = 0) : num_states(num_states), order(order) { }
		bool operator<(Index rhs) const
		{
			if (num_states < rhs.num_states)
			{
				return true;
			}
			return num_states == rhs.num_states && order < rhs.order;
		}
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & num_states;
			ar & order;
		}
	};
	typedef std::map<Index, model_ptr_t> model_map_t;



	/**
	Data
	*/

	/** The models. */
	model_map_t models;






	/**
	Methods
	*/

	DnaHmmOrderNumStateMap()
	{
	}

	/** Does the map have a particular model? */
	bool contains_model(unsigned num_states, unsigned order)
	{
		const Index index(num_states, order);
		model_map_t::iterator i = models.find(index);
		return models.end() != i;
	}

	/** Get the model for a particular order. */
	model_t & get_model(unsigned num_states, unsigned order)
	{
		return *(get_model_ptr(num_states, order));
	}

	/** Get the model for a particular order. */
	void insert_model(unsigned num_states, unsigned order, model_ptr_t model)
	{
		models[Index(num_states, order)] = model;
	}

	/** Get the model for a particular order. */
	model_ptr_t & get_model_ptr(unsigned num_states, unsigned order)
	{
		const Index index(num_states, order);
		model_map_t::iterator i = models.find(index);
		if (models.end() == i)
		{
			throw BIO_MAKE_STRING("Could not find model with " << num_states << " states of order " << order);
		}
		return i->second;
	}

	/** Train all the HMMs on the given sequences. */
	void train_all(const SequenceCollection & seq_list)
	{
		std::for_each(
			models.begin(),
			models.end(),
			ModelTrainer(seq_list));
	}

	void
	gen_sequence_from_random_hmm(seq_t & seq, size_t seq_length) const
	{
		if (0 == models.size())
		{
			throw std::logic_error( "Have not got any HMM DNA models" );
		}

		//choose an HMM to generate the sequence
		size_t index = get_uniform_index(models.size());
		model_map_t::const_iterator hmm = models.begin();
		while (0 != index)
		{
			++hmm;
			--index;
		}

		seq.clear();
		hmm->second->append_random_sequence(seq, seq_length);
	}





	/**
	Classes
	*/

	/**
	Trains a model given training data.
	*/
	struct ModelTrainer
	{
		const SequenceCollection & training_data;

		ModelTrainer(const SequenceCollection & training_data)
			: training_data(training_data)
		{
		}

		void operator()(model_map_t::value_type model)
		{
			model.second->train(training_data);
		}
	};



	void init_singleton();



	/**
	Serialisation
	*/
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & models;
    }

};




BIO_NS_END

#endif //BIO_HMM_DNA_H_
