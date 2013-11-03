#ifndef BIOBASE_HMM_SEQ_GEN_H_
#define BIOBASE_HMM_SEQ_GEN_H_

#include "bio/defs.h"
#include "bio/random.h"
#include "bio/hidden_markov_model.h"



BIO_NS_START

template <class Alphabet>
struct HmmSequenceGenerator
{
	typedef HiddenMarkovModel<Alphabet> hmm_t;

	HmmSequenceGenerator(const hmm_t * hmm)
		: hmm(hmm)
	{
		//choose a random initial state
		double rnd_value = get_uniform_01();
		for (state_idx = 0; rnd_value > hmm->states[state_idx].initial_prob; ++state_idx) {
			assert(state_idx < hmm->states.size());
			rnd_value -= hmm->states[state_idx].initial_prob;
		}
		assert(0 <= rnd_value);
	}

	Alphabet operator()()
	{
		//choose a random emission
		double rnd_value = get_uniform_01();
		size_t e;
		for (e = 0; rnd_value > hmm->states[state_idx].emission_probs[e]; ++e) {
			assert(e < hmm->states[state_idx].emission_probs.size());
			rnd_value -= hmm->states[state_idx].emission_probs[e];
		}
		assert(0 <= rnd_value);
		Alphabet generated_symbol = AlphabetTraits<Alphabet>::get_symbol(e);

		//choose a random transition
		rnd_value = get_uniform_01();
		size_t s;
		for (s = 0; rnd_value > hmm->states[state_idx].transition_probs[s]; ++s) {
			assert(s < hmm->states.size());
			rnd_value -= hmm->states[state_idx].transition_probs[s];
		}
		assert(0 <= rnd_value);
		state_idx = s;

		return generated_symbol;
	}

	const hmm_t * hmm;
	size_t state_idx;
};

BIO_NS_END




#endif //BIOBASE_HMM_SEQ_GEN_H_
