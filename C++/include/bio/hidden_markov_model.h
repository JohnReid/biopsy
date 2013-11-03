#ifndef BIO_HIDDEN_MARKOV_MODEL_H_
#define BIO_HIDDEN_MARKOV_MODEL_H_

#include "bio/defs.h"
#include "bio/hmm_alphabet.h"
#include "bio/random.h"

#include <boost/io/ios_state.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/access.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <algorithm>
#include <numeric>



BIO_NS_START

typedef std::vector<prob_t> prob_vector_t;
typedef std::vector<prob_vector_t> prob_matrix_t;



/**
Functor to test whether value is finite.
*/
struct NotFinite
{
	bool operator()(prob_t value) { return ! BIO_FINITE(value); }
};


/**
MarkovState represents one state in a hidden markov model.
*/
template <class Alphabet>
struct MarkovState {

	//
	// typedefs
	//
	typedef Alphabet alphabet_t;



	//
	// data
	//
	prob_t initial_prob; //the initial probability of this state
	prob_vector_t emission_probs; //the probability of emitting a symbol from the alphabet
	prob_vector_t transition_probs; //the probability of moving to another state



	//
	// methods
	//
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & initial_prob;
        ar & emission_probs;
        ar & transition_probs;
    }

	bool operator==(const MarkovState<Alphabet> & rhs) const {
		BIO_FPC_NS::close_at_tolerance<prob_t> test(1); //test at tolerance of 1%

		//check the emission probs are close
		prob_vector_t::const_iterator i = emission_probs.begin();
		prob_vector_t::const_iterator j = rhs.emission_probs.begin();
		for ( ; emission_probs.end() != i; ++i, ++j) {
			if (! test(*i, *j)) {
				return false;
			}
		}

		//check the transition probs are close
		i = transition_probs.begin();
		j = rhs.transition_probs.begin();
		for ( ; transition_probs.end() != i; ++i, ++j) {
			if (! test(*i, *j)) {
				return false;
			}
		}

		return test(initial_prob, rhs.initial_prob);
	}

	MarkovState () { }

	template <class EmIt, class TrIt>
	MarkovState<Alphabet> (
		prob_t initial_prob,
		EmIt emission_prob_begin,
		TrIt trans_begin,
		TrIt trans_end)
		: initial_prob(initial_prob)
	{
		std::copy(
			emission_prob_begin,
			emission_prob_begin + AlphabetTraits<alphabet_t>::get_size(),
			inserter(emission_probs, emission_probs.end()));

		std::copy(
			trans_begin,
			trans_end,
			inserter(transition_probs, transition_probs.end()));

		assert(is_consistent());
	}

	template <class A>
	prob_t
	get_emission_prob(A a) const
	{
		const unsigned idx = AlphabetTraits<alphabet_t>::get_index(a);
		if (idx > AlphabetTraits<alphabet_t>::get_size())
		{
			throw std::logic_error( "Index out of range" );
		}
		return emission_probs[idx];
	}

	//have default tolerance of 1 %
	bool is_consistent(prob_t tolerance = 0.01) const
	{
		if (NotFinite()(initial_prob)) {
			return false;
		}

		if (emission_probs.end() != std::find_if(emission_probs.begin(), emission_probs.end(), NotFinite())) {
			return false;
		}

		if (transition_probs.end() != std::find_if(transition_probs.begin(), transition_probs.end(), NotFinite())) {
			return false;
		}

		if (emission_probs.size() != AlphabetTraits<alphabet_t>::get_size()) {
			return false;
		}

		const prob_t sum = std::accumulate(emission_probs.begin(), emission_probs.end(), 0.0);
		if (1.0 - tolerance > sum || sum > 1.0 + tolerance) {
			return false;
		}
		return true;
	}
};
template <class OStr, class Alphabet>
OStr &
operator<<(OStr & os, const MarkovState<Alphabet> & state)
{
	boost::io::ios_precision_saver ips(os);
	os.precision(3);

	os << "Initial prob: " << state.initial_prob << std::endl;
	os << "Emission probabilities: ";
	for (prob_vector_t::const_iterator i = state.emission_probs.begin(); state.emission_probs.end() != i; ++i) {
		os << AlphabetTraits<Alphabet>::get_symbol(i - state.emission_probs.begin()) << "=" << *i << " ";
	}
	os << std::endl;
	os << "Transition probabilities: ";
	for (prob_vector_t::const_iterator i = state.transition_probs.begin(); state.transition_probs.end() != i; ++i) {
		os << (i - state.transition_probs.begin()) << "=" << *i << " ";
	}
	return os;
}





/**
HiddenMarkovModel
*/
template <class Alphabet>
struct HiddenMarkovModel
{
	//
	// typedefs
	//
	typedef Alphabet alphabet_t;
	typedef MarkovState<Alphabet> state_t;
	typedef std::vector<state_t> state_vector_t;
	typedef boost::shared_ptr<HiddenMarkovModel <Alphabet> > ptr_t;


	//
	// data
	//
	state_vector_t states; //the states


	//
	// methods
	//
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & states;
    }

	/** Default constructor. */
	HiddenMarkovModel<Alphabet> () { }

	/** Creates a HMM with the given number of states and randomised emission and transition probabilities. */
	HiddenMarkovModel<Alphabet>(size_t num_states)
	{
		prob_vector_t initial_probs;
		gen_prob_dist(num_states, inserter(initial_probs, initial_probs.begin()));

		for (size_t s = 0; s < num_states; ++s)
		{
			prob_vector_t emission_probs;
			gen_prob_dist(AlphabetTraits<Alphabet>::get_size(), inserter(emission_probs, emission_probs.begin()));

			prob_vector_t transition_probs;
			gen_prob_dist(num_states, inserter(transition_probs, transition_probs.begin()));

			states.push_back(
				state_t(
					initial_probs[s],
					emission_probs.begin(),
					transition_probs.begin(),
					transition_probs.end()));
		}

		assert(is_consistent());
	}

	bool operator==(const HiddenMarkovModel<Alphabet> & rhs) const { return states == rhs.states; }

	bool is_consistent(prob_t tolerance = 0.01) const
	{
		prob_t sum_initial_probs = 0.0;

		//first check each state is internally consistent and has the right number of transition probabilities
		for (typename state_vector_t::const_iterator i = states.begin();
			states.end() != i;
			++i)
		{
			if (! i->is_consistent()) {
				return false;
			}
			if (i->transition_probs.size() != states.size()) {
				return false;
			}
			sum_initial_probs += i->initial_prob;
		}

		if (1.0 - tolerance > sum_initial_probs || sum_initial_probs > 1.0 + tolerance) {
			return false;
		}

		return true;
	}
};
template <class Alphabet>
std::ostream &
operator<<(std::ostream & os, const HiddenMarkovModel<Alphabet> & hmm)
{
	typedef HiddenMarkovModel<Alphabet> hmm_t;
	for (typename hmm_t::state_vector_t::const_iterator i = hmm.states.begin();
		hmm.states.end() != i;
		++i)
	{
		os
			<< "State: " << i - hmm.states.begin() << std::endl
			<< *i
			<< std::endl;
	}
	return os;
}


BIO_NS_END

#endif //BIO_HIDDEN_MARKOV_MODEL_H_
