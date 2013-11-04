#
# Copyright John Reid 2008,2009
#

"""
Code to implement a MEME-like single gap motif finding algorithm.
"""


import hmm, hmm.pssm, logging, time, math, numpy, os, cPickle
from hmm.pssm import single_gap
from gapped_pssms.sequence import numpy_to_seq, seq_to_numpy
from gapped_pssms.initialisations import K_mer_distance, DistanceFilter, yield_k_mers
from gapped_pssms.mosaic import create_mosaic_model
from gapped_pssms.background import k_mer_log_likelihoods, forward_backward_log_likelihoods
from gapped_pssms.sites import sites_in_sequences, yield_sites_in_sequence, remove_sites
from gapped_pssms.output import make_output_dir
from gapped_pssms.pssm_score import calculate_first_order_entropy_score, calculate_information_content_score, geometric_mean
from itertools import ifilter, imap


def nucleo_dist_from_mer(mer, pseudo_count, gap_index=None):
    # initialise our distribution to 'N's
    mer_len = len(mer)
    nucleo_dist = numpy.ones((mer_len+1, 4)) * pseudo_count

    # add K_mer to distribution
    for i, c in enumerate(mer):
        #from IPython.Debugger import Pdb; Pdb().set_trace()
        if None != gap_index and i >= gap_index: # skip the gap if necessary
            index = i+1
        else:
            index = i
        if 4 == c:
            nucleo_dist[index] += .25
        else:
            nucleo_dist[index,c] += 1

    # normalise
    for dist in nucleo_dist:
        dist /= dist.sum()

    return nucleo_dist

def evaluate_gap_position(
  L_mer,
  gap_index,
  sequences,
  bg_L_mer_scores,
  pssm_scores,
  comp_pssm_scores,
  options
):
    "Evaluate a k-mer with a gap in a particular position."
    gap_dist = nucleo_dist_from_mer(L_mer, options.pseudo_count_for_L_mer_scoring, gap_index=gap_index)
    gap_pssm = hmm.calculate_log_scores(gap_dist)
    gap_comp_pssm = hmm.calculate_complementary_scores(gap_pssm)
    gap_pssm_scores = hmm.max_scores_in_sequences(gap_pssm, sequences, bg_L_mer_scores)
    gap_comp_pssm_scores = hmm.max_scores_in_sequences(gap_comp_pssm, sequences, bg_L_mer_scores)
    max_scores_per_pssm = numpy.array(
        [
            pssm_scores,
            comp_pssm_scores,
            gap_pssm_scores,
            gap_comp_pssm_scores
        ]
    )
    best_scores = max_scores_per_pssm.max(axis=0)
    best_scores[best_scores<0.] = 0. # when we didn't find sites, ignore
    score = best_scores.sum() / len(L_mer) / len(best_scores)
    logging.debug(
        'Evaluated: %s; gap: %d; score: %f',
        numpy_to_seq(L_mer),
        gap_index,
        score
    )
    return score

def evaluate_L_mer(L_mer, sequences, bg_L_mer_scores, gap_positions, options):
    dist = nucleo_dist_from_mer(
      L_mer,
      options.pseudo_count_for_L_mer_scoring,
      gap_index=None
    )
    pssm = hmm.calculate_log_scores(dist)
    comp_pssm = hmm.calculate_complementary_scores(pssm)
    pssm_scores = hmm.max_scores_in_sequences(pssm, sequences, bg_L_mer_scores)
    comp_pssm_scores = hmm.max_scores_in_sequences(comp_pssm, sequences, bg_L_mer_scores)
    return max(
      (
        evaluate_gap_position(
          L_mer,
          gap_index,
          sequences,
          bg_L_mer_scores,
          pssm_scores,
          comp_pssm_scores,
          options
        ),
        gap_index
      )
      for gap_index in gap_positions
    )


def generate_seeds(
  sequences,
  preprocessed_sequences,
  options
):
    """
    Generate a list of candidate L-mers and score them to find the best seed L-mer and gap position.
    """
    # if we have been given a background model filename and it exists then load it.
    if None != options.bg_model_filename and os.path.exists(options.bg_model_filename):
        logging.info("Loading supplied background model from %s", options.bg_model_filename)
        bg_model = cPickle.load(open(options.bg_model_filename))
        converted_seqs = [bg_model.converter.to_order_n(s) for s in sequences]
    else:
        logging.info("Learning new background model")
        bg_model, converted_seqs = learn_bg_model(
          sequences,
          num_mosaics=options.bg_model_num_mosaics,
          order=options.bg_model_order
        )
        if options.bg_model_filename:
            logging.info("Saving background model to %s", options.bg_model_filename)
            cPickle.dump(bg_model, open(options.bg_model_filename, 'w'))

    if options.force_seed:
        logging.info('Forcing seed to be: %s', options.force_seed)
        L_mers = [(seq_to_numpy(options.force_seed), len(sequences), len(sequences))]
    else:
        # Calculate log likelihood of L-mers under background model.
        bg_L_mer_scores = calculate_k_mer_scores(bg_model, converted_seqs, options.L)

        # Find best candidate L-mers
        distance = K_mer_distance(allowed_shifts=options.allowed_shifts, shift_cost=options.shift_cost)
        L_mer_seeds = list()
        gap_end_offset = options.L/5 + 1
        start = time.time()
        num_L_mers_to_find = 3 * options.max_L_mers_to_evaluate
        logging.info('Finding best %d candidate %d-mers to seed HMM emissions', num_L_mers_to_find, options.L)
        L_mers = hmm.top_mers_by_sequence_membership(
          preprocessed_sequences,
          k=options.L,
          n=num_L_mers_to_find
        )
        logging.info('Finding top %d %d-mers took %f seconds', len(L_mers), options.L, time.time()-start)

    if options.force_gap:
        logging.info('Forcing gap at position: %d', options.force_gap)
        L_mer_seeds = [
          (numpy_to_seq(L_mer), L_mer_count, L_mer_num_seqs, options.force_gap, 0.0)
          for L_mer, L_mer_count, L_mer_num_seqs in L_mers
        ]
    else:
        # Evaluate L-mers
        if -1 == options.seed_filter_distance:
            min_distance = options.L / 4 + 1
        else:
            min_distance = options.seed_filter_distance
        logging.info('Positioning gaps up to %d bases from end of K-mers', gap_end_offset)
        gap_positions = range(gap_end_offset, options.L+1-gap_end_offset)
        logging.info('Filtering K-mers that are not %d away from previously evaluated.', min_distance)
        logging.info('Evaluating up to %d L-mers.', options.max_L_mers_to_evaluate)
        L_mer_filter = DistanceFilter(distance, min_distance=min_distance)
        discarded = 0
        evaluated = 0
        for L_mer, L_mer_count, L_mer_num_seqs in L_mers:
            if not L_mer_filter(L_mer) or 4 in L_mer:
                logging.debug('Discarding: %s; count: %d; # sequences: %d', numpy_to_seq(L_mer), L_mer_count, L_mer_num_seqs)
                discarded += 1
            else:
                score, gap_index = evaluate_L_mer(L_mer, sequences, bg_L_mer_scores, gap_positions, options)
                evaluated += 1
                logging.info(
                    'Evaluated (%3d/%d): %s; gap: %d; count: %d; # sequences: %d; score: %f',
                    evaluated, options.max_L_mers_to_evaluate, numpy_to_seq(L_mer), gap_index, L_mer_count, L_mer_num_seqs, score
                )
                L_mer_seeds.append((numpy_to_seq(L_mer), L_mer_count, L_mer_num_seqs, gap_index, score))
                if len(L_mer_seeds) == options.max_L_mers_to_evaluate:
                    break
        L_mer_seeds.sort(key=lambda x: -x[4]) # sort by score, highest first
        logging.info('Discarded %d L-mers using edit distance', discarded)
        logging.info('Evaluated %d L-mers: scores range from %f to %f', evaluated, L_mer_seeds[-1][4], L_mer_seeds[0][4])
    return L_mer_seeds, bg_model




class SingleGapAlgorithm(object):
    """
    MEME-like single gap motif finding algorithm.
    """

    def __init__(
      self,
      options,
    ):
        self.options = options
        "The parameters for the algorithm."

    def __call__(self, sequences, bg_model=None):
        """
        Run the motif finding algorithm.
        """
        logging.info('Looking for at least %d PSSMs', self.options.num_pssms)

        if self.options.max_L_mers_to_evaluate < self.options.num_pssms:
            raise ValueError('Cannot find any more PSSMs than L-mers are evaluated.')

        num_bases = sum(len(s) for s in sequences)
        logging.info('Running single gap algorithm on %d sequences with %d bases', len(sequences), num_bases)

        preprocessed_sequences = [hmm.preprocess_sequence(s) for s in sequences]
        start = time.clock()
        seeds, bg_model = generate_seeds(
          sequences,
          preprocessed_sequences,
          self.options
        )
        logging.info('Generating %d seeds took %.1f seconds', len(seeds), time.clock() - start)

        # try the best few seeds
        start = time.clock()
        num_to_examine = self.options.num_seeds_to_examine and self.options.num_seeds_to_examine or 2 * self.options.num_pssms
        logging.info('Examining %d/%d seeds', num_to_examine, len(seeds))
        p_one_per_seq = len(preprocessed_sequences) / float(num_bases) # expect one per sequence
        self.p_binding_site = p_one_per_seq * self.options.p_binding_site_scale
        logging.info(
            'HMM p(binding site) parameter estimated as %f (1 site/seq) and adjusted to %f by scaling parameter',
            p_one_per_seq,
            self.p_binding_site,
        )
        results = list(
          (seed, self.try_seed(seed, len(sequences), num_bases, preprocessed_sequences))
          for seed, i in zip(seeds, xrange(num_to_examine))
        )

        # keep only those results that succeeded
        results = filter(lambda x: x[1], results)
        logging.info('Got PSSMs for %d seeds in %.1f seconds', len(results), time.clock() - start)

        # define a function that scores results
        def score_result(result):
            seed, result = result
            L_mer, count, num_seqs, gap, score = seed
            model, builder, num_sites, num_seqs_with_site = result
            emissions, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)
            return geometric_mean((
              calculate_first_order_entropy_score(emissions),
              calculate_information_content_score(emissions),
              num_seqs_with_site / float(len(sequences))
            ))

        # sort by scores
        scored_results = [(score_result(result), result) for result in results]
        scored_results.sort(reverse=True)

        # remove those PSSMs that do not score highly enough
        logging.info('Removing PSSMs with low scores.')
        while len(scored_results) > self.options.num_pssms and scored_results[-1][0] < self.options.pssm_score_threshold:
            scored_results.pop()
        logging.info('%d PSSMs scored highly enough', len(scored_results))

        # examine results
        for i, (score, (seed, result)) in enumerate(scored_results):
            logging.info('************** PSSM %d **************', i)
            L_mer, count, num_seqs, gap, seed_score = seed
            model, builder, num_sites, num_seqs_with_site = result
            logging.info(
              'Seed %s with gap at %d had %d hits in %d/%d sequences',
              L_mer,
              gap,
              count,
              num_seqs,
              len(sequences)
            )
            logging.info('Seed score: %f', seed_score)
            image_file = os.path.join(self.options.output_dir, '%s-%03d' % (self.options.tag, i))
            pssm_def_file = os.path.join(self.options.output_dir, '%s-%03d.pssm' % (self.options.tag, i))
            logging.info(
              'HMM found %d sites. %d/%d sequences have at least one site',
              num_sites,
              num_seqs_with_site,
              len(sequences)
            )
            self.examine_model(model, builder, sequences, image_file, pssm_def_file)
            logging.info('Score: %g', score)
            emissions, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)

        return [r[1] for r in results]


    def try_seed(self, seed, num_seqs, num_bases, preprocessed_sequences):
        """
        Try the given L-mer as a starting point for Baum-Welch.
        """
        best_L_mer, best_L_mer_count, best_L_mer_num_seqs, best_gap_index, best_score = seed
        logging.info(
          'Trying seed: %s; gap: %2d; # hits: %5d; # sequences: %4d; score: %f',
          best_L_mer,
          best_gap_index,
          best_L_mer_count,
          best_L_mer_num_seqs,
          best_score
        )
        #from IPython.Debugger import Pdb; Pdb().set_trace()

        #
        # Build our model
        #
        start, builder = self._make_builder(best_gap_index, len(best_L_mer))
        model = self._model_for_L_mer(
          best_L_mer,
          best_gap_index,
          self.p_binding_site
        )

        #
        # How many sites does the initialised model find?
        #
#    start = time.time()
#    num_sites, num_seqs_with_site = sites_in_sequences(model, (0,), preprocessed_sequences)
#    logging.info('Took %f seconds to find sites in %d bases', time.time()-start, num_bases)
#    logging.info('Initial model found %d sites. %d/%d sequences have at least one site', num_sites, num_seqs_with_site, num_seqs)
#    # if we cannot find sites in more than a given fraction of sequences we discard this seed
#    fraction_found = num_seqs_with_site/float(num_seqs)
#    if num_seqs_with_site/float(num_seqs) < self.options.initialised_model_threshold:
#      logging.info('Not enough sites (%f). Skipping this seed', fraction_found)
#      return None

        #
        # Train our model
        #
        LLs = list()
        def _callback(LL):
            "A callback for every step of Baum-Welch."
            LLs.append(LL)
            emissions, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)
            entropy_per_base = hmm.pssm.entropy(emissions, gap_probs) / gap_probs.sum()
            logging.debug('LL: %f; Entropy/base: %f', LL, entropy_per_base)
            return bool(entropy_per_base < self.options.entropy_per_base_threshold)
        tolerance = self.options.bw_tolerance_per_seq*num_seqs
        prior = hmm.ModelPrior(model.N, model.M)
        A = numpy.zeros_like(prior.A)
        A[0] = model.A[0] * num_bases * 5
        if self.options.force_no_gap:
            # we've been asked to not use the gap so place a large prior on the transitions that avoid it.
            gi = builder.gap_index
            A[gi,gi+2] = 1e9
            rgi = builder.K + gi + 2
            A[rgi,rgi-2] = 1e9
        prior.A = A
        logging.debug('Running Baum-Welch with tolerance: %f', tolerance)
        start = time.time()
        LL, num_iterations = model.baum_welch(
          preprocessed_sequences,
          tolerance=tolerance,
          prior=prior,
          callback=_callback,
        )
        logging.info('Baum-Welch achieved LL %f in %d iterations (%f secs)', LL, num_iterations, time.time()-start)
        emissions, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)
        entropy_per_base = hmm.pssm.entropy(emissions, gap_probs) / gap_probs.sum()
        if entropy_per_base >= self.options.entropy_per_base_threshold:
            logging.info('Baum-Welch halted early as entropy/base = %f > %f threshold', entropy_per_base, self.entropy_per_base_threshold)
            return None

        #
        # How many sites does the trained model find?
        #
        sites = [list(yield_sites_in_sequence(model, s, (0,))) for s in preprocessed_sequences]
        num_sites, num_seqs_with_site = sum(len(s) for s in sites), sum(len(s) > 0 for s in sites)
        logging.info('Trained model found %d sites. %d/%d sequences have at least one site', num_sites, num_seqs_with_site, len(preprocessed_sequences))

        #logging.info('Blanking sites')
        #for i, (pp_seq, site) in enumerate(zip(preprocessed_sequences, sites)):
        #  seq = pp_seq.as_numpy()
        #  # convert to correct type
        #  seq_int = numpy.asarray(seq, dtype=int)
        #  #from IPython.Debugger import Pdb; Pdb().set_trace()
        #  remove_sites(site, seq_int, unknown_output=model.M)
        #  preprocessed_sequences[i] = hmm.preprocess_sequence(seq_int)

        return model, builder, num_sites, num_seqs_with_site


    def examine_model(self, model, builder, sequences, image_file=None, pssm_def_file=None):
        """
        Log some info about the model.
        """
        #
        # How many sites does it find after training?
        #
        emissions, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)
        logging.info('Entropy/base        : %f', hmm.pssm.entropy(emissions, gap_probs) / gap_probs.sum())
        logging.info('Information content : %f', hmm.pssm.information_content(emissions))
        if None != pssm_def_file:
            output_pssm_definition(open(pssm_def_file, 'w'), emissions, gap_probs)
        if None != image_file:
            import hmm.pssm.logo as logo
            image = logo.pssm_as_image(emissions, transparencies=gap_probs)
            png_file = '%s.png' % image_file
            logging.info('Saving PSSM to %s', png_file)
            image.save(png_file, "PNG")
            eps_file = '%s.eps' % image_file
            logging.info('Saving PSSM to %s', eps_file)
            image.save(eps_file, "EPS")


    def _make_builder(self, gap_index, mer_len):
        # calculate where to put the gap in the longer PSSM
        start = (self.options.K - (mer_len+1))/2
        return start, single_gap.SingleGappedPssmBuilder(
          K=self.options.K,
          gap_index=gap_index+start,
          markov_order=0,
          M=4
        )


    def _model_for_L_mer(self, L_mer, gap_index, p_binding_site):
        """
        Create a model initialised by this K-mer.
        """
        # get the start position of the K-mer and a builder to make the model
        mer_len = len(L_mer)
        start, builder = self._make_builder(gap_index, mer_len)

        # get the emission distribution
        nucleo_dist = nucleo_dist_from_mer(
          seq_to_numpy(L_mer),
          self.options.pseudo_count_for_model_initialisation,
          gap_index=gap_index
        )
        emissions = numpy.ones((self.options.K,4))/4.
        emissions[start:start+mer_len+1] = nucleo_dist

        # build the model
        pssm, in_states, out_states = builder.create(
          p_gap=.5,
          emissions=emissions
        )
        model = hmm.as_model(
          single_gap.add_to_simple_background_model(
            model=pssm,
            in_states=in_states,
            out_states=out_states,
            p_binding_site=p_binding_site
          )
        )
        #print model.A
        #from IPython.Debugger import Pdb; Pdb().set_trace();
        return model



def output_pssm_definition(f, emissions, gap_probs):
    K = len(emissions)
    assert len(gap_probs) == K
    print >> f, "# Output of output_pssm_definition()"
    print >> f, "# in %s" % __file__
    print >> f, "MODEL:%d,%d" % (K, emissions.shape[1])
    print >> f, "INITIAL:0;1"
    for i in xrange(1, K):
        gp = gap_probs[i]
        print >> f, "TRANSITION:%d;%d;%g" % (i-1, i, gp)
        if gp < 1.:
            print >> f, "TRANSITION:%d;%d;%g" % (i-1, i+1, 1.-gp)
    for i, e in enumerate(emissions):
        print >> f, "EMISSIONS:%d;%s" % (i, ','.join(imap(str, e)))

def uniform_bg_model():
    bg_model = hmm.as_model(
      create_mosaic_model(
        num_mosaics=1,
        p_transition=.1,
        alphabet_size=4,
        order=0,
        dirichlet_prior_strength=None
      )
    )
    return bg_model


def learn_bg_model(
  sequences,
  num_mosaics=4,
  order=3,
  tolerance_per_base=7e-5
):
    """
    @return: (bg_model, converted_sequences)
    """
    bg_model = hmm.as_model(
      create_mosaic_model(
        num_mosaics=4,
        p_transition=.1,
        alphabet_size=4,
        order=3,
        dirichlet_prior_strength=.3
      )
    )
    converted_seqs = [bg_model.converter.to_order_n(s) for s in sequences]
    def _callback(LL):
        logging.debug('Background model LL: %f', LL)
        return True
    tolerance=tolerance_per_base*sum(len(s) for s in sequences)
    logging.info('Learning background model with %d mosaics of order %d, tolerance=%f', num_mosaics, order, tolerance)
    start = time.time()
    LL, iterations = bg_model.baum_welch(converted_seqs, tolerance=tolerance, callback=_callback)
    logging.info('Achieved LL=%f after %d iterations and %f seconds', LL, iterations, time.time() - start)
    return bg_model, converted_seqs



def calculate_k_mer_scores(bg_model, converted_seqs, K):
    result = list()
    for seq in converted_seqs:
        LL, alpha, beta, c = bg_model.forward_backward(seq)
        k_mer_LLs = k_mer_log_likelihoods(K=K, LL=LL, alpha=alpha, beta=beta, c=c)
        result.append(k_mer_LLs)
    return result

def add_algorithm_options(option_parser):
    """
    Add options to an option parser that control the parameters of the algorithm.
    """
    option_parser.add_option(
      "-o",
      "--output",
      dest="output_dir",
      default='.',
      help="The directory to write output to."
    )
    option_parser.add_option(
      "-b",
      "--bg-model-filename",
      dest="bg_model_filename",
      help="The filename of the pickled background model."
    )
    option_parser.add_option(
      "-t",
      "--tag",
      dest="tag",
      default='single-gap',
      help="String to identify this run of the single gap algorithm. Used as part of output filenames."
    )
    option_parser.add_option(
      "-K",
      dest="K",
      type='int',
      default=16,
      help="Length of the PSSMs to look for."
    )
    option_parser.add_option(
      "-L",
      dest="L",
      type='int',
      default=8,
      help="Length of the L-mers to seed the HMM with."
    )
    option_parser.add_option(
      "--seed-filter-distance",
      dest="seed_filter_distance",
      type='int',
      default=-1,
      help="Minimum Hamming distance between candidate seeds."
    )
    option_parser.add_option(
      "-N",
      dest="num_pssms",
      type='int',
      default=1,
      help="Find at least this number of PSSMs."
    )
    option_parser.add_option(
      "--force-seed",
      dest="force_seed",
      type='string',
      default=None,
      help="Force algorithm to use given L-mer as seed."
    )
    option_parser.add_option(
      "--force-gap",
      dest="force_gap",
      type='int',
      default=None,
      help="Force algorithm to use given gap position."
    )
    option_parser.add_option(
      "--no-gap",
      dest="force_no_gap",
      action='store_true',
      default=False,
      help="Force algorithm to not use any gap."
    )
    option_parser.add_option(
      "--num-seeds-to-examine",
      dest="num_seeds_to_examine",
      type='int',
      default=0,
      help="Number of seeds to examine. If 0, examines twice as many seeds as -N option."
    )
    option_parser.add_option(
      "--max-L-mers-to-evaluate",
      dest="max_L_mers_to_evaluate",
      type='int',
      default=8,
      help="Maximum number of L-mers to evaluate."
    )
    option_parser.add_option(
      "--pssm-score-threshold",
      dest="pssm_score_threshold",
      type='float',
      default=.5,
      help="PSSMs will be rejected unless they score at least this highly."
    )
    option_parser.add_option(
      "--bg-model-num-mosaics",
      dest="bg_model_num_mosaics",
      type='int',
      default=4,
      help="Number of mosaics (states) in HMM to model background."
    )
    option_parser.add_option(
      "--bg-model-order",
      dest="bg_model_order",
      type='int',
      default=3,
      help="Markov order of HMM to model background."
    )
    option_parser.add_option(
      "--allowed-shifts",
      dest="allowed_shifts",
      type='int',
      default=2,
      help="Number of shifts allowed when filtering L-mers."
    )
    option_parser.add_option(
      "--shift-cost",
      dest="shift_cost",
      type='int',
      default=1,
      help="Cost of each shift when filtering L-mers."
    )
    option_parser.add_option(
      "--pseudo-count-for-model-initialisation",
      dest="pseudo_count_for_model_initialisation",
      type='float',
      default=.05,
      help="The strength of the pseudo count used when we create a PSSM from a K-mer to initialise our HMM."
    )
    option_parser.add_option(
      "--pseudo-count-for-L-mer-scoring",
      dest="pseudo_count_for_L_mer_scoring",
      type='float',
      default=.01,
      help="The strength of the pseudo count used when we create a PSSM to score the L-mer seeds."
    )
    option_parser.add_option(
      "--entropy-per-base-threshold",
      dest="entropy_per_base_threshold",
      type='float',
      default=1.25,
      help="Baum-Welch will halt if entropy/base is higher than this."
    )
    option_parser.add_option(
      "--initialised-model-threshold",
      dest="initialised_model_threshold",
      type='float',
      default=.01,
      help="If a model initialised from a seed does not find sites in this fraction of sequences it is discarded before Baum-Welch training."
    )
    option_parser.add_option(
      "--bw_tolerance_per_seq",
      dest="bw_tolerance_per_seq",
      type='float',
      default=4e-5,
      help="Stopping criteria for change in log likelihood when running Baum-Welch."
    )
    option_parser.add_option(
      "--p-binding-site-scale",
      dest="p_binding_site_scale",
      type='float',
      default=.05,
      help="Scales the estimate of p(binding site) in the HMMs. Can be used to tune for more specific or more vague motifs."
    )
