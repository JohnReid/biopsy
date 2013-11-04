#
# Copyright John Reid 2010
#


"""
Code to compare scoring methods in the BiFa module.
"""


import oreganno_dataset as dataset, biopsy as B, logging, pylab as P, random as R, numpy as N

from infpy import roc


def make_sequence_vec(strings):
    "@return: Return the collection of strings into a biopsy.SequenceVec."
    result = B.SequenceVec()
    result.extend(strings)
    return result




class BiFaClassifier(object):
    "Classifies test cases using BiFa algorithm."

    def __init__(self, min_threshold=1e-2):
        self.min_threshold = min_threshold

    def score(self, test_case):
        "Score the test case and return the hits."
        record, seq, pssms = test_case
        hits = B.score_pssms_on_sequence(make_sequence_vec(pssms), str(seq.seq), threshold=self.min_threshold)
        return hits

    def __call__(self, test_case):
        "@return: The maximum threshold at which this test case would be classified as positive."
        hits = self.score(test_case)
        if not hits:
            return self.min_threshold
        else:
            return max(hit.p_binding for hit in hits)




def thresholds_for_classifier(test_cases, classifier):
    "Classify the test cases and return the positive and negative thresholds."
    positive_thresholds = []
    negative_thresholds = []
    for test_case, truth in test_cases:
        threshold = classifier(test_case)
        if truth:
            positive_thresholds.append(threshold)
        else:
            negative_thresholds.append(threshold)

    positive_thresholds.sort()
    negative_thresholds.sort()

    return positive_thresholds, negative_thresholds




def plot_param_setting_rocs(test_cases, use_cumulative_dists, use_p_value, use_score, label, color, linestyle):
    "Test a particular setting of the biopsy module's PSSM parameters."
    B.PssmParameters.use_cumulative_dists = use_cumulative_dists
    B.PssmParameters.use_p_value = use_p_value
    B.PssmParameters.use_score = use_score
    start_time = time.time()
    positive_thresholds, negative_thresholds = thresholds_for_classifier(test_cases, BiFaClassifier())
    elapsed = time.time() - start_time
    logging.info(
        'Testing parameter settings: %40s : %10.3f : %10.1f secs',
        label,
        sum(positive_thresholds) + sum(negative_thresholds),
        elapsed
    )
    rocs = roc.rocs_from_thresholds(positive_thresholds, negative_thresholds, num_points=100)
    roc.plot_roc_points(rocs, label=label, color=color, ls=linestyle)
    return positive_thresholds, negative_thresholds





logging.basicConfig(level=logging.INFO)

R.seed(1)

param_settings = (
    (True,  False, True,  'Cumulative Bayesian',        'blue', '-.'),
    (False, False, True,  'Non-cumulative Bayesian',    'blue', '-'),
    (False, False, False, 'BiFa',                       'cyan', '-'),
    (True,  True,  True,  'Cumulative p-values',        'red',  '-.'),
    (False, True,  True,  'Non-cumulative p-values',    'red',  '-'),
)


try: test_cases
except NameError:
    test_cases = list(dataset.generate_test_cases(ignore_chip_data=False))
logging.info('Have %d test cases', len(test_cases))


#P.close('all')
P.figure(figsize=(16,12))
for param_setting in param_settings:
    plot_param_setting_rocs(test_cases, *param_setting)

roc.plot_random_classifier()
P.xlim(0, 0.2)
P.ylim(0, 0.4)
P.xlim(0, 1)
P.ylim(0, 1)
P.legend(loc='lower right')
P.show()
