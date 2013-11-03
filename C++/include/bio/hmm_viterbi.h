#ifndef BIO_VITERBI_H_
#define BIO_VITERBI_H_

#include "bio/defs.h"
#include "bio/hidden_markov_model.h"
USING_BIO_NS

#include <vector>
#include <limits>
using namespace std;

#include <cmath>

#ifdef min
# undef min
#endif
#ifdef max
# undef max
#endif


BIO_NS_START

/** See "A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition." Lawrence R. Rabiner
 use_logs specifies whether we multiply probabilities or add their logarithms
 */
template <bool use_logs>
struct ViterbiAlgorithm
{
	typedef vector<size_t> index_vector_t;
	typedef vector<index_vector_t> psi_matrix_t;
	psi_matrix_t psi;
	prob_matrix_t delta; 

	/**
	HMM is a hidden markov model
	ObsIt is an iterator over the observed sequence
	InsIt is a front_insert_iterator for the most likely sequence of states
	*/
	template <
		class HMM,
		class ObsIt,
		class InsIt>
	void viterbi(
		const HMM & hmm,
		size_t num_obs,
		ObsIt o_begin,
		ObsIt o_end,
		InsIt state_insert_it)
	{
		//delta stores the best score (probability) along a single path which accounts for the first
		//so many observations and ends in a given state
		delta.resize(num_obs);

		//psi stores the states that maximised delta for given t and j
		psi.resize(num_obs);

		//do we have anything to do?
		if (0 == num_obs) {
			//no
			return;
		}

		const size_t num_states = hmm.states.size();

		//initialisation
		size_t t = 0;
		size_t i = 0;
		psi[t].resize(num_states);
		delta[t].resize(num_states);
		for (typename HMM::state_vector_t::const_iterator s = hmm.states.begin(); hmm.states.end() != s; ++s, ++i) {
			delta[t][i] =
				use_logs
					? std::log(s->initial_prob) + std::log(s->get_emission_prob(*o_begin))
					: s->initial_prob * s->get_emission_prob(*o_begin);
		}
		++o_begin;
		++t;

		//recursion - for each observed emission
		for ( ; o_end != o_begin; ++o_begin) {

			//put another row in our matrices
			psi[t].resize(num_states);
			delta[t].resize(num_states);

			//for each state
			size_t j = 0; //the state index
			for (typename HMM::state_vector_t::const_iterator s = hmm.states.begin(); hmm.states.end() != s; ++s, ++j) {

				//find the maximum - we update these 2 variables as we find the maximum over all the states
				prob_t max_prob = -std::numeric_limits<prob_t>::max();
				size_t max_index = 0;
				
				i = 0;
				for (typename HMM::state_vector_t::const_iterator s2 = hmm.states.begin(); hmm.states.end() != s2; ++s2, ++i) {
					prob_t new_prob = numeric_limits<prob_t>::min();
					new_prob =
						use_logs 
							? delta[t-1][i] + std::log(s2->transition_probs[j])
							: delta[t-1][i] * s2->transition_probs[j];
					if (new_prob > max_prob) {
						max_prob = new_prob;
						max_index = i;
					}
				}

				psi[t][j] = max_index;
				delta[t][j] =
					use_logs
						? max_prob + std::log(s->get_emission_prob(*o_begin))
						: max_prob * s->get_emission_prob(*o_begin);
			}

			++t;
		}

		assert(t == num_obs);

		//termination
		size_t q_star = 0;
		for (i = 1; i < num_states; ++i) {
			if (delta[t-1][i] > delta[t-1][q_star]) {
				q_star = i;
			}
		}
		*state_insert_it++ = q_star;

		//path (state sequence) backtracking
		for ( ; t != 0; --t) {
			q_star = psi[t-1][q_star];
			*state_insert_it++ = q_star;
		}
	}
};


BIO_NS_END

#endif //BIO_VITERBI_H_

