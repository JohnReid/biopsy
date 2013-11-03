#ifndef BIO_FORWARD_BACKWARD_H_
#define BIO_FORWARD_BACKWARD_H_

#include "bio/hidden_markov_model.h"

#include <deque>
#include <algorithm>

#include <cmath>

BIO_NS_START


/**
See "A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition." Lawrence R. Rabiner

The use_scaling parameter specifies if we use the scaling implementation designed to avoid underflows

*/
template <bool use_scaling>
struct ForwardBackwardAlgorithm
{
	prob_matrix_t alpha; //forward matrix
	prob_matrix_t beta; //backward matrix
	prob_vector_t c; //scaling coefficients

	/**
	HMM is a hidden markov model
	ObsIt is an iterator over the observed sequence
	*/
	template <
		class HMM,
		class ObsIt
    >
	ForwardBackwardAlgorithm< use_scaling > &
	forward(
		const HMM & hmm,
		ObsIt o_begin,
		ObsIt o_end)
	{
		alpha.clear();
		c.clear();

		const size_t num_states = hmm.states.size();
		size_t i, j;

		//do we have anything to do?
		if (o_begin != o_end) {

			//initialisation
			size_t t = 0;
			alpha.push_back(prob_vector_t());
			alpha[t].resize(num_states);
			i = 0;
			for (i = 0; num_states != i; ++i) {
				alpha[t][i] = hmm.states[i].initial_prob * hmm.states[i].get_emission_prob(*o_begin);
			}
			if (use_scaling) {
				c.push_back(1.0);
			}
			++o_begin;
			++t;

			//induction - i.e. for each observation
			for ( ; o_end != o_begin; ++o_begin) {

				alpha.push_back(prob_vector_t());
				alpha[t].resize(num_states);
				assert(alpha.size() - 1 == t);

				for (j = 0; num_states != j; ++j) {

					prob_t sum = 0.0;

					for (i = 0; num_states != i; ++i) {
						sum += alpha[t-1][i] * hmm.states[i].transition_probs[j];
					}

					alpha[t][j] = sum * hmm.states[j].get_emission_prob(*o_begin);

				}

				if (use_scaling)
				{

					//scale if necessary for floating point precision
					const prob_t sum = std::accumulate(alpha[t].begin(), alpha[t].end(), 0.0);

					if (0.0 != sum && sum < 0.001) //do we want to scale?
					{
						
						//scale
						c.push_back(1.0 / sum);
						assert(BIO_FINITE(c[t]));
						for (i = 0; i < num_states; ++i) {
							alpha[t][i] *= c[t];
						}

					} else {
						//no scaling
						c.push_back(1.0);
					}
				}

				++t;

			} //for each observation

		} //if o_begin != o_end

		if (use_scaling) {
			assert(c.size() == alpha.size());
		}

		return *this;
	}

	/**
	HMM is a hidden markov model
	ObsRIt is a reverse iterator over the observed sequence
	*/
	template <
		class HMM,
		class ObsRIt>
	void
	backward(
		const HMM & hmm,
		ObsRIt o_rbegin,
		ObsRIt o_rend)
	{
		const size_t num_obs = alpha.size();
		beta.resize(num_obs);

		const size_t num_states = hmm.states.size();
		size_t i, j;

		//do we have anything to do?
		if (o_rbegin != o_rend) {

			size_t t = num_obs - 1;
			beta[t].resize(num_states);
			std::fill(beta[t].begin(), beta[t].end(), 1.0); //beta[T-1] is an array of 1's
			--t;

			//induction - i.e. for each observation - we don't use the first (last in reverse order) observation
			while (o_rbegin + 1 != o_rend)
			{

				beta[t].resize(num_states);

				for (i = 0; i != num_states; ++i)
				{

					prob_t sum = 0.0;
					for (j = 0; j != num_states; ++j)
					{
						const prob_t p_transition = hmm.states[i].transition_probs[j];
						const prob_t p_emission = hmm.states[i].get_emission_prob(*o_rbegin);
						const prob_t p_state_at_t_plus_1 = beta[t+1][j];

						BOOST_ASSERT(BIO_FINITE(p_transition));
						BOOST_ASSERT(BIO_FINITE(p_emission));
						BOOST_ASSERT(BIO_FINITE(p_state_at_t_plus_1));

						sum +=
							p_transition
							* p_emission
							* p_state_at_t_plus_1;

					} //for each state, j
					assert(BIO_FINITE(sum));

					beta[t][i] = sum;

				} //for each state, i

				if (use_scaling)
				{

					for (size_t i = 0; i < num_states; ++i)
					{
						
						assert(BIO_FINITE(c[t]));
						beta[t][i] *= c[t]; //use scaling factor
						if (! BIO_FINITE(beta[t][i]))
						{
							throw std::logic_error( "Overflow" );
						}

					} //for each state, i

				} //if (use_scaling)

				--t;
				++o_rbegin;

			} //for each observation

		} //if o_rbegin != o_rend
	}

	/** get the probability of the sequence used to run the forward algorithm. */
	prob_t
	get_probability()
	{
		if (0 == alpha.size()) {
			//empty sequences have probability 1
			return 1.0;
		}

		//termination
		const prob_t unscaled_result = std::accumulate(alpha.rbegin()->begin(), alpha.rbegin()->end(), 0.0);
		if (use_scaling) {
			const prob_t scale = std::accumulate(c.begin(), c.end(), 1.0, std::multiplies<prob_t>());
			assert(0.0 != scale);
			return unscaled_result / scale;
		} else {
			return unscaled_result;
		}
	}

	/** get the log probability of the sequence used to run the forward algorithm. */
	prob_t
	get_log_probability()
	{
		if (0 == alpha.size())
		{
			//empty sequences have probability 1
			return std::log(1.0);
		}

		//termination
		const prob_t unscaled_result = std::accumulate(alpha.rbegin()->begin(), alpha.rbegin()->end(), 0.0);

		if (use_scaling)
		{
			prob_t result = std::log(unscaled_result);
			for (prob_vector_t::const_iterator i = c.begin(); c.end() != i; ++i)
			{
				assert(BIO_FINITE(*i));

				result -= std::log(*i);
			}
			return result;
		}
		else
		{
			return std::log(unscaled_result);
		}
	}
};

/** Get the log likelihood of a sequence given the hmm. */
template <class HMM, class SeqIt>
prob_t
get_log_likelihood(const HMM & hmm, SeqIt seq_begin, SeqIt seq_end)
{
	return ForwardBackwardAlgorithm<true>()
		.forward(
			hmm,
			seq_begin,
			seq_end).get_log_probability();
}

/** Get the average likelihood of a character in the sequence given the hmm. */
template <class HMM, class SeqIt>
prob_t
get_per_char_likelihood(const HMM & hmm, SeqIt seq_begin, SeqIt seq_end)
{
	const prob_t log_prob = get_log_likelihood(hmm, seq_begin, seq_end);

	const unsigned size = unsigned(seq_end - seq_begin);

	return (0 == size) ? 1.0 : std::exp(log_prob / size);
}

/** Get the accumulated log likelihood of every sequence given the hmm. */
template <class HMM, class SeqListIt>
prob_t
get_multiple_log_likelihood(const HMM & hmm, SeqListIt seq_list_begin, SeqListIt seq_list_end)
{
	prob_t sum = 0;
	for ( ; seq_list_end != seq_list_begin; ++seq_list_begin)
	{
		sum += get_log_likelihood(hmm, seq_list_begin->begin(), seq_list_begin->end());
	}
	return sum;
}

BIO_NS_END

#endif //BIO_FORWARD_BACKWARD_H_
