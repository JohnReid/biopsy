#ifndef BIO_BAUM_WELCH_H_
#define BIO_BAUM_WELCH_H_

#include "bio/hmm_forward_backward.h"

#include <algorithm>
#include <limits>
#include <functional>

#ifdef min
# undef min
#endif
#ifdef max
# undef max
#endif


BIO_NS_START



template< bool use_scaling >
struct BaumWelchAlgorithm : public ForwardBackwardAlgorithm< use_scaling >
{
	/** Numerator of the estimate of the transition probablities. */
	prob_matrix_t xi_sum_t;

	/** Numerator of the estimate of the emission probablities. */
	prob_matrix_t gamma_sum_obs_v;

	/** Denominator of the estimate of the emission and transition probablities. */
	prob_vector_t gamma_sum;

	/** The transition probabilities to update to. */
	prob_matrix_t transition_probs;

	/** The initial estimates to update to. */
	prob_vector_t pi;

	/**
	HMM is a hidden markov model
	ObsIt is an iterator over the observed sequence
	*/
	template <
		class HMM,
		class ObsIt,
		class ObsRIt>
	void
	run(
		HMM & hmm,
		ObsIt o_begin,
		ObsIt o_end,
		ObsRIt o_rbegin,
		ObsRIt o_rend)
	{
		if (o_begin != o_end)
		{
			//separate calculation and updating so that multiple algorithm can just call first
			calculate_estimates(
				hmm,
				o_begin,
				o_end,
				o_rbegin,
				o_rend);

			update(hmm);
		}
	}

	/**
	HMM is a hidden markov model
	ObsIt is an iterator over the observed sequence
	*/
	template <
		class HMM,
		class ObsIt,
		class ObsRIt>
	void
	calculate_estimates(
		HMM & hmm,
		ObsIt o_begin,
		ObsIt o_end,
		ObsRIt o_rbegin,
		ObsRIt o_rend)
	{
		const size_t num_states = hmm.states.size();
		const size_t alphabet_size = AlphabetTraits<typename HMM::alphabet_t>::get_size();

		ForwardBackwardAlgorithm< use_scaling >::template forward(
			hmm,
			o_begin,
			o_end);

		ForwardBackwardAlgorithm< use_scaling >::template backward(
			hmm,
			o_rbegin,
			o_rend);

		const size_t num_obs = this->alpha.size();
		assert(this->beta.size() == num_obs);

		//stores the probability of a pair of consecutive states
		prob_matrix_t xi_t; //for current t
		xi_t.resize(num_states);
		xi_sum_t.resize(num_states);
		for (size_t i = 0; i < num_states; ++i) //resize and initialise matrices
		{
			xi_t[i].resize(num_states);
			xi_sum_t[i].resize(num_states);
			std::fill(xi_sum_t[i].begin(), xi_sum_t[i].end(), 0.0);
		}

		//stores the sum of a row of xi
		prob_vector_t gamma; //gamma for current t
		gamma.resize(num_states);
		gamma_sum.resize(num_states);
		std::fill(gamma_sum.begin(), gamma_sum.end(), 0.0);
		gamma_sum_obs_v.resize(num_states);
		for (size_t i = 0; i < num_states; ++i)
		{
			gamma_sum_obs_v[i].resize(alphabet_size);
			std::fill(gamma_sum_obs_v[i].begin(), gamma_sum_obs_v[i].end(), 0.0);
		}

		//stores the newly calculated transition probabilities
		transition_probs.clear();
		transition_probs.resize(num_states);

		//for each observation
		ObsIt o = o_begin;
		size_t t;
		for (t = 0; t != num_obs; ++o)
		{

			prob_t xi_sum = 0.0;
			std::fill(gamma.begin(), gamma.end(), 0.0);

			//for each state, i
			for (size_t i = 0; i < num_states; ++i)
			{

				//for each state, j
				for (size_t j = 0; j < num_states; ++j)
				{
					const prob_t p_state_i_at_t = this->alpha[t][i];
					const prob_t p_transition_to_j = hmm.states[i].transition_probs[j];
					const prob_t p_observe_o = hmm.states[i].get_emission_prob(*o);
					const prob_t p_state_j_at_t_plus_1 = this->beta[t][j];

					BOOST_ASSERT(BIO_FINITE(p_state_i_at_t));
					BOOST_ASSERT(BIO_FINITE(p_transition_to_j));
					BOOST_ASSERT(BIO_FINITE(p_observe_o));
					BOOST_ASSERT(BIO_FINITE(p_state_j_at_t_plus_1));

					xi_t[i][j] =
						p_state_i_at_t
						* p_transition_to_j
						* p_observe_o
						* p_state_j_at_t_plus_1;

					xi_sum += xi_t[i][j];
					BOOST_ASSERT(BIO_FINITE(xi_sum));

				} //for each state, j

			} //for each state, i

			//deal with zero denominator
			if (0.0 == xi_sum)
			{
				//used to throw here but it seems legitimate that a sequence might have 0
				//likelihood, so...
				xi_sum = std::numeric_limits<prob_t>::min();

				//if xi_sum is 0 so will all gamma[i]
				//throw std::logic_error( "Probabilities vanished" );
			
			}

			//divide each element by the sum and update xi_sum_t
			for (size_t i = 0; i < num_states; ++i)
			{
				for (size_t j = 0; j < num_states; ++j)
				{

					xi_t[i][j] /= xi_sum;
					if (! BIO_FINITE(xi_t[i][j]))
					{
						throw std::logic_error( "Overflow" );
					}

					xi_sum_t[i][j] += xi_t[i][j];
					gamma[i] += xi_t[i][j];

				} //for each state, j

				//update gamma_sum and gamma_sum_obs_v
				gamma_sum[i] += gamma[i];

				const size_t index_of_last_symbol = AlphabetTraits<typename HMM::alphabet_t>::get_index(*o);
				gamma_sum_obs_v[i][index_of_last_symbol] += gamma[i];

			} //for each state, i

			//if it is the first observation we can estimate the initial probabilities
			if (0 == t)
			{
				const prob_t initial_prob_sum = std::accumulate(gamma.begin(), gamma.end(), 0.0);
				if (0.0 == initial_prob_sum)
				{
					//if sum(gamma[i]) is 0 then distribute the initial states evenly
					std::fill(gamma.begin(), gamma.end(), 1.0 / gamma.size());
				}
				else if (0.999 > initial_prob_sum || initial_prob_sum > 1.001)
				{
					//normalise
					std::transform(
						gamma.begin(),
						gamma.end(),
						gamma.begin(),
						std::bind1st(std::multiplies<prob_t>(), 1.0 / initial_prob_sum));
				}

				pi.clear();
				for (size_t i = 0; i < num_states; ++i)
				{
					pi.push_back(gamma[i]);
				} //for each state, i

			}

			++t;

			//if we have got to T-1 save the transition probabilities
			if (t == num_obs - 1)
			{
				for (size_t i = 0; i < num_states; ++i)
				{
					if (gamma_sum[i] == 0.0)
					{
						//no reason to update emission or transition probs as we have not been in state i
						//we leave the i'th vector in 'a' empty as a marker of this fact.
					}
					else
					{
						//transition probs
						for (size_t j = 0; j < num_states; ++j)
						{
							const prob_t new_transition_prob = xi_sum_t[i][j] / gamma_sum[i];
							if (! BIO_FINITE(new_transition_prob))
							{
								throw std::logic_error( "Overflow" );
							}
							transition_probs[i].push_back(new_transition_prob);

						} //for each state, j

					}
				}
			}

		} //for each time, t
	}

	template <class HMM>
	void update(HMM & hmm)
	{
		const size_t num_states = hmm.states.size();
		const size_t alphabet_size = AlphabetTraits<typename HMM::alphabet_t>::get_size();

		//for each state
		for (size_t i = 0; i != num_states; ++i)
		{
			//update initial probs
			hmm.states[i].initial_prob = pi[i];

			//only update transition probs if we have something to update with
			if (! transition_probs[i].empty())
			{
				//can update transition probs here as well
				for (size_t j = 0; j < num_states; ++j)
				{
					hmm.states[i].transition_probs[j] = transition_probs[i][j];
				}
			}

			//only update emission probs if we visited the state
			if (gamma_sum[i] != 0.0)
			{
				for (size_t k = 0; k < alphabet_size; ++k)
				{
					const prob_t new_emission_prob = gamma_sum_obs_v[i][k] / gamma_sum[i];
					if (! BIO_FINITE(new_emission_prob))
					{
						throw std::logic_error( "Overflow" );
					}
					hmm.states[i].emission_probs[k] = new_emission_prob;

				} //for each alphabet symbol, k
			}
		}
	}
};



/** Executes the baum welch algorithm for the multiple sequence case. */
template <bool use_scaling>
struct BaumWelchMultipleAlgorithm
{
	typedef BaumWelchAlgorithm<use_scaling> alg_t;
	typedef std::vector<alg_t> alg_vec;

	alg_vec algs;


	template <
		class HMM,
		class SeqIt>
	void
	run(
		HMM & hmm,
		SeqIt seq_begin,
		SeqIt seq_end)
	{
		calculate_estimates(hmm, seq_begin, seq_end);

		update(hmm);
	}


	template <
		class HMM,
		class SeqIt>
	void
	calculate_estimates(
		HMM & hmm,
		SeqIt seq_begin,
		SeqIt seq_end)
	{
		//construct a vector of single sequence baum welch algorithms
		for (SeqIt s = seq_begin;
			seq_end != s;
			++s)
		{
			algs.push_back(alg_t());
			algs.rbegin()->calculate_estimates(
				hmm,
				s->begin(),
				s->end(),
				s->rbegin(),
				s->rend());
		}
	}



	template <class HMM>
	void
	update(HMM & hmm)
	{
		const size_t num_states = hmm.states.size();
		const size_t alphabet_size = AlphabetTraits<typename HMM::alphabet_t>::get_size();

		//update the transition and emission probabilities
		for (size_t i = 0; i < num_states; ++i)
		{
			//the initial probs
			prob_t initial_prob = 0;
			unsigned num_samples = 0;
			for (typename alg_vec::const_iterator a = algs.begin();
				algs.end() != a;
				++a)
			{
				if (! a->pi.empty()) //can be empty if empty test sequence
				{
					initial_prob += a->pi[i];
					++num_samples;
				}
			}
			if (0 != num_samples)
			{
				hmm.states[i].initial_prob = initial_prob / num_samples;
			}


			//transition probs
			for (size_t j = 0; j < num_states; ++j)
			{
				prob_t numerator = 0.0;
				unsigned num_samples = 0;

				for (typename alg_vec::const_iterator a = algs.begin();
					algs.end() != a;
					++a)
				{
					//were we in this state in this sequence?
					if (a->transition_probs[i].empty())
					{
						//no - so there is no contribution from this sequence to this state
						continue;
					}

					numerator += a->transition_probs[i][j];
					++num_samples;
				}

				if (0.0 != num_samples)
				{
					hmm.states[i].transition_probs[j] = numerator / num_samples;
					if (! BIO_FINITE(hmm.states[i].transition_probs[j]))
					{
						throw std::logic_error( "Overflow" );
					}
				}
			}

			//emission probs
			for (size_t l = 0; l < alphabet_size; ++l)
			{
				prob_t numerator = 0.0;
				prob_t denominator = 0.0;

				for (typename alg_vec::const_iterator a = algs.begin();
					algs.end() != a;
					++a)
				{
					//were we in this state in this sequence?
					if (0 == a->gamma_sum[i])
					{
						//no - so there is no contribution from this sequence to this state
						continue;
					}

					numerator += a->gamma_sum_obs_v[i][l];
					denominator += a->gamma_sum[i];
				}

				if (0.0 != denominator)
				{
					hmm.states[i].emission_probs[l] = numerator / denominator;
					if (! BIO_FINITE(hmm.states[i].emission_probs[l]))
					{
						throw std::logic_error( "Overflow" );
					}
				}
			} //for each alphabet symbol, l
		}
	}
};

template <class HMM, class SeqIt, class SeqRIt>
void
baum_welch_single(HMM & hmm, SeqIt seq_begin, SeqIt seq_end, SeqRIt seq_rbegin, SeqRIt seq_rend)
{
	BaumWelchAlgorithm<true>().run(
		hmm,
		seq_begin,
		seq_end,
		seq_rbegin,
		seq_rend);
}

template <class HMM, class SeqListIt>
void
baum_welch_multiple(HMM & hmm, SeqListIt seq_list_begin, SeqListIt seq_list_end)
{
	BaumWelchMultipleAlgorithm<true>().run(
		hmm,
		seq_list_begin,
		seq_list_end);
}



BIO_NS_END

#endif //BIO_BAUM_WELCH_H_
