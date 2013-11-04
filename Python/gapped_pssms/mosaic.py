#
# Copyright John Reid 2008
#

"""
Code to create mosaic models.
"""

import hmm, hmm.pssm, numpy

def create_mosaic_model(
  num_mosaics,
  p_transition,
  alphabet_size,
  order,
  dirichlet_prior_strength=None
):
    """
    Create a mosaic model.

    Each mosaic has an independent parameter that specifies the probability of transitioning to any other
    given mosaic (this effectively ties these transition probabilities together).

    @arg num_mosaics: The number of mosaics.
    @arg p_transition: The probability of leaving a mosaic.
    @arg alphabet_size: The size of the output alphabet.
    @arg order: The Markov order of the output.
    @arg dirichlet_prior_strength: The strength of the uniform prior on the emission probabilities. If None
    then a uniform distribution is used.
    """

    builder = hmm.pssm.ModelBuilder(order, alphabet_size=alphabet_size)
    model = builder.new_model_by_states()

    for n in xrange(num_mosaics):
        # add the state
        if None == dirichlet_prior_strength:
            emission_dist = numpy.ones(model.M) / alphabet_size
        else:
            emission_dist = hmm.dirichlet_draw(numpy.ones(model.M)*dirichlet_prior_strength)
        state = builder.add_fully_parameterised_state(
          model,
          pi=1./num_mosaics,
          emission_dist=emission_dist
        )

    if 1 == num_mosaics:
        model.states[0].add_successor(model.states[0], model.add_parameter(1.))
    else:
        for s1 in model.states:
            transition_param = model.add_parameter(p_transition)
            no_transition_param = model.add_parameter(1.0 - p_transition)
            for s2 in model.states:
                s1.add_successor(s2, s1 == s2 and no_transition_param or transition_param)

    return model

def evaluate_mosaics(max_mosaics=6, max_order=3):
    """
    Evaluate different mosaic models on chip-chip fragments.
    """
    from gapped_pssms import data
    sequences = data.training_test_sequences()
    mosaic_sizes = range(1,max_mosaics+1)
    orders = range(max_order+1)
    preprocessed_sequences = [
      (
        [hmm.preprocess_sequence(s) for s in training],
        [hmm.preprocess_sequence(s) for s in test]
      )
      for training, test in sequences
    ]
    result = list()
    for order in orders:
        converter = hmm.MarkovOrderConverter(alphabet_size=4, order=order)
        order_n_seqs = [
          (
            [converter.to_order_n(s) for s in training],
            [converter.to_order_n(s) for s in test]
          )
          for training, test in sequences
        ]
        for num_mosaics in mosaic_sizes:
            LL = 0.
            for training_seqs, test_seqs in order_n_seqs:
                model = hmm.as_model(
                  create_mosaic_model(
                    num_mosaics=num_mosaics,
                    p_transition=0.,
                    alphabet_size=4,
                    order=order,
                    dirichlet_prior_strength=10.
                  )
                )
                model.baum_welch(training_seqs)
                LL += sum(model.LL(s) for s in test_seqs)
            logging.info('Order: %d; # mosaics: %d; LL: %f', order, num_mosaics, LL)
            result.append((order, num_mosaics, LL))
    return result



if '__main__' == __name__:
    from gapped_pssms.data import fasta_file_for_fragment
    from gapped_pssms.sequence import convert_fasta_sequences
    import logging
    from pylab import *

    logging.basicConfig(level=logging.INFO)

    max_mosaics = 15
    max_order = 5
    me = evaluate_mosaics(max_mosaics=max_mosaics, max_order=max_order)
    for i in range(max_order):
        e=me[i*max_mosaics:(i+1)*max_mosaics]
        plot([x[1] for x in e], [x[2] for x in e])
    xlabel('# mosaics')
    ylabel('LL')
    title('Evaluation of mosaic models of various Markov orders')
    savefig('mosaic-evaluation.png', format='PNG')
    raise

    # load our sequences
    sequences = convert_fasta_sequences(fasta_file_for_fragment('T00671'))

    # build our model
    model_by_states = create_mosaic_model(
      num_mosaics=1,
      p_transition=0.,
      alphabet_size=4,
      order=2,
      dirichlet_prior_strength=10.
    )
    model = hmm.as_model(model_by_states)
    print model.B

    # convert our sequences to the correct order
    sequences_order_n = [model.converter.to_order_n(s) for s in sequences]

    #from IPython.Debugger import Pdb; Pdb().set_trace()
    def callback(LL):
        logging.info('LL: %f', LL)
    model.baum_welch(sequences_order_n, callback=callback)
    print model.B
