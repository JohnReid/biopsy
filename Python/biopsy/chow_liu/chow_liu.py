#
# Copyright John Reid 2008
#

"""
Code to implement feature dependencies using the Chow-Liu algorithm.
Chow, C. K. & C. N. Liu (1968),
"Approximating discrete probability distributions with dependence trees",
IEEE Transactions on Information Theory IT-14 (3): 462-467.

Each I{datum} in a set of I{data} is represented as a list. Each entry in the list represents the value
of a I{feature} for that datum. Missing values are represented by None.
"""

import numpy as N, random, os, sys, fpconst
from cookbook.dicts import DictOf
from cookbook.pre_post_conditions import *


def is_finite_condition(x, *in_args):
    "Tests if first argument is finite. Used as post-condition for other methods."
    assert fpconst.isFinite(x)


def normalise_array(array):
    "@return: array/array.sum()"
    return array / array.sum()



def make_zero_where_arg_is_zero(result, arg):
    """
    @return: result amended such that it is 0 where arg == 0
    """
    result[N.where(0 == arg)] = 0
    return result



@postcondition(is_finite_condition)
def mutual_information(p_x_y):
    "@return: The mutual information between x and y (using a base 2 logarithm)."
    p_x_p_y = N.outer(normalise_array(p_x_y.sum(axis=1)), normalise_array(p_x_y.sum(axis=0)))
    return make_zero_where_arg_is_zero(
      p_x_y * (N.log2(p_x_y) - N.log2(p_x_p_y)),
      p_x_y
    ).sum()



def max_mutual_information(num_values):
    """
    @return: The maximum mutual information possible for the given number of values.

    TODO: Fix this! Just guess work for moment.
    """
    return N.log2(num_values)






def count_joint_features(data, num_features, num_values, pseudo_count=.35):
    """
    Counts the co-occurence of values of pairs of features in the data.

    @arg data: Each datum in data is a list of integers. Each element in each list is the value
    of the corresponding feature or None for missing data.
    @return: A map from tuples of features (f1, f2) to arrays of counts (c1, c2).
    """

    def initial_counts(key):
        f1, f2 = key
        if f1 == f2:
            return N.zeros(num_values) + pseudo_count
        else:
            return N.zeros((num_values, num_values)) + 4. * pseudo_count

    result = DictOf(initial_counts, take_key_as_arg=True)

    for datum in data:
        if num_features != len(datum):
            raise RuntimeError('Datum has wrong number (%d) of features (should be %d)' % (len(datum), num_features))
        for f1, v1 in enumerate(datum):
            if None != v1:
                assert v1 < num_values
                # update marginal prob of f1 feature
                result[f1,f1][v1] += 1.
                for f2, v2 in enumerate(datum[:f1]):
                    if None != v2:
                        assert v2 < num_values
                        result[f1,f2][v1,v2] += 1.

    return result


def counts_to_joint_distribution(counts):
    """
    Turns counts into joint probabilities.
    """
    return dict((features, normalise_array(c)) for features, c in counts.iteritems())


def expand_joint_dists(joint_dists, num_features):
    "Expands joint distributions to include views that represent dists with axes swapped."
    for f1 in xrange(num_features):
        for f2 in xrange(f1):
            if (f1, f2) in joint_dists:
                joint_dists[f2,f1] = joint_dists[f1,f2].swapaxes(0,1)


def mutual_informations_from_joint_dists(joint_dists):
    """
    @return: A list of tuples: ((f1, f2), joint, mi) where f1 and f2 are features, joint is their
    joint distribution and mi is their mutual information.
    """
    return [
      ((f1, f2), joint, mutual_information(joint))
      for (f1, f2), joint in joint_dists.iteritems()
      if f1 != f2
    ]



def sort_mutual_informations(mutual_infos):
    """
    Sort a list of mutual information tuples produced by L{mutual_informations_from_joint_dists} in-place.
    """
    def MI(mi):
        return mi[2]
    mutual_infos.sort(key=MI)


def strongest_dependency(mutual_infos):
    """
    @return: the dependency with the highest mutual information.
    """
    def MI(x):
        return x[2]
    return max(mutual_infos, key=MI)



def build_weighted_graph(mutual_infos, num_features, num_values):
    """
    @arg mutual_infos: The output from L{mutual_informations_from_joint_dists}.
    @arg num_features: The number of features in the data.
    @arg num_values: The number of values each feature can take.
    @return: B{g, features, weights} where
      - B{g} is a boost.graph.Graph
      - B{features} is a map from vertices to features
      - B{weights} is a map from edges to floats that can be used in a minimum spanning
      tree algorithm to find the Chow-Liu tree.
    """
    # build a graph and label each vertex with its feature
    import boost.graph as bgl
    g = bgl.Graph()
    vertices = [g.add_vertex() for f in xrange(num_features)]
    features = g.add_vertex_property(name='features', type='integer')
    for i, v in enumerate(vertices):
        features[v] = i

    # add the edges and label each edge with its mutual information
    weights = g.add_edge_property(type='float')
    mutual_info_map = g.add_edge_property(name='mutual infos', type='float')
    for (f1, f2), joint, mi in mutual_infos:
        #print f1, f2, mi
        e = g.add_edge(vertices[f1], vertices[f2])
        mutual_info_map[e] = mi
        weights[e] = max_mutual_information(num_values) - mi
        assert weights[e] >= 0.
    return g, features, weights


def compute_chow_liu_tree(g, weights, features, min_mutual_info=0.):
    """
    @arg g: The graph of all mutual informations.
    @arg weights: Map from edges to weights for minimum spanning tree algorithm.
    @arg features: Map from vertices to features.
    @return: A map from features to their predecessors in the minimum spanning tree.
    """
    import boost.graph as bgl
    # calculate minimal spanning tree on graph
    predecessor_map = g.add_vertex_property(type='vertex')
    mst_edges = bgl.prim_minimum_spanning_tree(g, predecessor_map=predecessor_map, weight_map=weights)

    # calucate the predecessor of each feature
    feature_predecessors = dict()
    mutual_infos = g.edge_properties['mutual infos']
    for v in g.vertices:
        if v != predecessor_map[v]:
            e = g.edge(predecessor_map[v], v)
            if mutual_infos[e] > min_mutual_info:
                feature_predecessors[features[v]] = features[predecessor_map[v]]

    return feature_predecessors


def highlight_tree_edges(g, features, predecessor_map):
    """
    Change the style of the edges to dashed for those edges not in the Chow-Liu tree.

    @arg g: The graph of features and their mutual information weights.
    @arg features: Map from vertices to features.
    @arg predecessor_map: The map that defines the Chow-Liu tree.
    """
    style = g.add_edge_property(name='style', type='string')
    for e in g.edges:
        f1, f2 = features[g.source(e)], features[g.target(e)]
        if (f2 in predecessor_map and predecessor_map[f2] == f1) \
          or \
          (f1 in predecessor_map and predecessor_map[f1] == f2):
            style[e] = 'solid'
        else:
            style[e] = 'dashed'




def write_graph(g, features, weights, name, num_values):
    """
    Write a graph to graphviz file and convert it to SVG.
    """
    vertex_labels = g.add_vertex_property(name='label', type='string')
    for v in g.vertices:
        vertex_labels[v] = '%d' % features[v]
    edge_labels = g.add_edge_property(name='label', type='string')
    for e in g.edges:
        edge_labels[e] = '%.2g' % (max_mutual_information(num_values) - weights[e])
    dot_filename = '%s.dot' % name
    g.write_graphviz(dot_filename)
    os.system('dot %s.dot -T svg -o %s.svg' % (name, name))




def calculate_conditioners(g, feature_predecessors, num_features):
    """
    Calculates the conditioners (i.e. the feature that each feature is conditioned on) in the Chow-Liu tree
    defined by the predecessor map, feature_predecessors.
    """
    conditioners = [None] * num_features
    for target, source in feature_predecessors.iteritems():
        conditioners[target] = source
    return conditioners




def conditional_from_joint(joint_dist):
    """
    @arg joint_dist: joint_dist[f,c] = p(f,c)
    @return: cd where cd[f,c] = p(f|c)
    """
    return joint_dist/joint_dist.sum(axis=0)




def conditional_dist(feature, conditioner, joint_dists):
    """
    @return: m where m[f,c] = p(f|c)
    """
    assert feature != conditioner
    return conditional_from_joint(joint_dists[feature, conditioner])




def conditional_log_dist(feature, conditioner, joint_dists):
    """
    @return: None if conditioner is None, otherwise the log condition dist, m where m[f,c] = log p(f|c)
    """
    if None == conditioner:
        return None
    return N.log(conditional_dist(feature, conditioner, joint_dists))




def calculate_log_probs(conditioners, joint_dists):
    """
    Calculates the marginal log probabilities of each feature's values and also the conditional
    log probabilities for the predecessors given in the predecessor map.
    """
    log_marginals = [
      N.log(joint_dists[f,f])
      for f in xrange(len(conditioners))
    ]
    log_conditionals = [
      conditional_log_dist(feature, conditioner, joint_dists)
      for feature, conditioner in enumerate(conditioners)
    ]
    return log_marginals, log_conditionals



def calculate_marginal_log_likelihood(datum, log_marginals):
    """
    The log likelihood of the datum under the feature independent model.

    @return: The marginal likelihood of the datum given the log marginals of the feature values.
    """
    return sum(
      log_marginals[f][v]
      for f, v in enumerate(datum)
      if None != v
    )

def calculate_chow_liu_log_likelihood(datum, conditioners, log_marginals, log_conditionals):
    """
    @return: The likelihood of the datum given the Chow-Liu tree. I.e. using the conditional
    likelihoods.
    """
    result = 0.0
    for f, v in enumerate(datum):
        if None != v:
            conditioner = conditioners[f]
            if None == conditioner or None == datum[conditioner]:
                result += log_marginals[f][v]
            else:
                result += log_conditionals[f][v,datum[conditioner]]
    return result




def calculate_log_likelihood_ratio(datum, conditioners, log_marginals, log_conditionals):
    """
    @return: The log likelihood ratio between the Chow-Liu tree model and the independent feature model
    defined by the log marginals.
    """
    return calculate_chow_liu_log_likelihood(datum, conditioners, log_marginals, log_conditionals) \
      - calculate_marginal_log_likelihood(datum, log_marginals)




class DependencyAnalyser(object):
    """
    Analyses the dependencies in the data.
    """

    def __init__(self, data, num_features, num_values, pseudo_count=.35, min_mutual_info=0.):
        self.num_features = num_features
        self.num_values = num_values
        self.pseudo_count = pseudo_count
        self.min_mutual_info = min_mutual_info
        self.counts = count_joint_features(data, num_features, num_values, pseudo_count=pseudo_count)
        self.joint_dists = counts_to_joint_distribution(self.counts)
        self.mutual_infos = mutual_informations_from_joint_dists(self.joint_dists)
        self.g, self.features, self.weights = build_weighted_graph(self.mutual_infos, num_features, num_values)
        self.feature_predecessors = compute_chow_liu_tree(
          self.g,
          self.weights,
          self.features,
          min_mutual_info=self.min_mutual_info
        )
        self.conditioners = calculate_conditioners(self.g, self.feature_predecessors, num_features)
        expand_joint_dists(self.joint_dists, num_features)
        self.log_marginals, self.log_conditionals = calculate_log_probs(self.conditioners, self.joint_dists)

    def highlight_tree_edges(self):
        highlight_tree_edges(self.g, self.features, self.feature_predecessors)

    def write_graph(self, filename):
        write_graph(self.g, self.features, self.weights, filename, self.num_values)

    def strongest_dependency(self):
        return strongest_dependency(self.mutual_infos)

    def calculate_marginal_log_likelihood(self, datum):
        return calculate_marginal_log_likelihood(datum, self.log_marginals)

    def calculate_chow_liu_log_likelihood(self, datum):
        return calculate_chow_liu_log_likelihood(datum, self.conditioners, self.log_marginals, self.log_conditionals)

    def calculate_mixed_log_likelihood(self, datum, alpha):
        """
        @return: chow_liu_log_likelihood when alpha = 0., marginal_log_likelihood when alpha = 1. or
        a weighted average when alpha is in-between.
        """
        if 0. > alpha or 1. < alpha:
            raise RuntimeError('alpha must be between 0 and 1')
        return N.log(
          alpha * N.exp(self.calculate_marginal_log_likelihood(datum))
          + (1. - alpha) * N.exp(self.calculate_chow_liu_log_likelihood(datum))
        )

    def remove_edges_below(self, mutual_info_threshold):
        mutual_infos = self.g.edge_properties['mutual infos']
        to_remove = [e for e in self.g.edges if mutual_infos[e] < mutual_info_threshold]
        for e in to_remove:
            self.g.remove_edge(e)



def predictive_log_likelihood_ratio(training, test, num_features, num_values):
    """
    Calculate the log likelihood ratio of the test data for feature-dependent and independent models
    trained on the training data.
    """
    if not training or not test:
        return 0.
    dependencies = DependencyAnalyser(training, num_features, num_values, pseudo_count=.35, min_mutual_info=0.)

    return sum(
      calculate_log_likelihood_ratio(
        datum,
        dependencies.conditioners,
        dependencies.log_marginals,
        dependencies.log_conditionals
      )
      for datum in test
    )



def test_pseudo_counts(data_sets, num_values, pseudo_counts, min_mutual_info=0., num_folds=3):
    """
    Look at predictive marginal and conditional log likelihoods for different pseudo counts.
    """
    for pseudo_count in pseudo_counts:
        marginal_ll = 0.
        conditional_ll = 0.
        for data in data_sets:
            shuffled_data = list(data)
            random.shuffle(shuffled_data)
            for fold in xrange(num_folds):
                training = [datum for i, datum in enumerate(data) if fold != i % num_folds]
                test = [datum for i, datum in enumerate(data) if fold == i % num_folds]
                if not training or not test:
                    continue
                dependencies = DependencyAnalyser(
                  training,
                  len(training[0]),
                  num_values,
                  pseudo_count=pseudo_count,
                  min_mutual_info=min_mutual_info
                )
                marginal_ll += sum(
                  calculate_marginal_log_likelihood(datum, dependencies.log_marginals)
                  for datum in test
                )
                conditional_ll += sum(
                  calculate_chow_liu_log_likelihood(
                    datum,
                    dependencies.conditioners,
                    dependencies.log_marginals,
                    dependencies.log_conditionals
                  )
                  for datum in test
                )
        yield pseudo_count, marginal_ll, conditional_ll


def test_mixed_log_likelihoods(data_sets, num_values, alphas, pseudo_count=.3, min_mutual_info=0., num_folds=3):
    """
    Look at predictive marginal and conditional log likelihoods for different pseudo counts.
    """
    for alpha in alphas:
        LL = 0.
        for data in data_sets:
            shuffled_data = list(data)
            random.shuffle(shuffled_data)
            for fold in xrange(num_folds):
                training = [datum for i, datum in enumerate(data) if fold != i % num_folds]
                test = [datum for i, datum in enumerate(data) if fold == i % num_folds]
                if not training or not test:
                    continue
                dependencies = DependencyAnalyser(
                  training,
                  len(training[0]),
                  num_values,
                  pseudo_count=pseudo_count,
                  min_mutual_info=min_mutual_info
                )
                LL += sum(
                  dependencies.calculate_mixed_log_likelihood(datum, alpha)
                  for datum in test
                )
        yield alpha, LL



def cross_validate(
  data,
  num_features,
  num_values,
  num_folds=3
):
    """
    Build feature-dependent and independent models of cross-folds of the data and compare
    predictive log likelihood.
    """
    data = list(data)
    random.shuffle(data)
    return sum(
      predictive_log_likelihood_ratio(
        [datum for i, datum in enumerate(data) if fold != i % num_folds],
        [datum for i, datum in enumerate(data) if fold == i % num_folds],
        num_features,
        num_values
      )
      for fold in xrange(num_folds)
    )


if '__main__' == __name__:
    from read_transfac_sites import *
    import pylab

    # test sequences
    if True:
        seqs = [
          'acgt',
          'ttgt',
          'acaa',
          'ttaa'
        ]
        num_features = len(seqs[0])
        num_values = 4
        data = map(convert_seq, seqs)
        dependencies = DependencyAnalyser(data, num_features, num_values, pseudo_count=.0, min_mutual_info=.2)
        dependencies.highlight_tree_edges()
        dependencies.write_graph('seqs')
        test_seq = 'acaa'
        converted_test = convert_seq(test_seq)
        print '%s: marginal=%6.6f conditional=%6.6f' % (
          test_seq,
          dependencies.calculate_marginal_log_likelihood(converted_test),
          dependencies.calculate_chow_liu_log_likelihood(converted_test),
        )


    # test transfac sites
    if True:
        print 'Loading transfac sites'
        transfac_sites = load_transfac_sites()
        converted_sites = convert_transfac_sites(transfac_sites)
        print 'Loaded sites for %d matrices' % len(converted_sites)

        if False:
            print 'Testing different pseudo-counts'
            for pseudo_count, marginal_ll, conditional_ll in test_pseudo_counts(
              converted_sites.values(),
              4,
              [.001, .25, .3, .325, .35, .375, 1.],
              min_mutual_info=0.
            ):
                print '%4.3f %6.1f %6.1f' % (pseudo_count, marginal_ll, conditional_ll)


        if True:
            print 'Testing different weighted mixtures of conditional and marginal log likelihoods'
            alpha = N.arange(0, 1.05, .1)
            for alpha, LL in test_mixed_log_likelihoods(
              converted_sites.values(),
              4,
              alphas=alpha,
              pseudo_count=.3,
              min_mutual_info=0.
            ):
                print '%4.3f %6.1f' % (alpha, LL)


        if True:
            num_seqs = []
            LL_ratios = []
            best_mis = []
            for matrix, seqs in converted_sites.iteritems():
                try:
                    num_seq = len(seqs)
                    LL_ratio = cross_validate(seqs, len(seqs[0]), 4) / num_seq
                    dependencies = DependencyAnalyser(seqs, len(seqs[0]), 4)
                    print '%s %3d sequences: %6.6g' % (matrix, num_seq, LL_ratio)
                    LL_ratios.append(LL_ratio)
                    num_seqs.append(num_seq)
                    best_mis.append(dependencies.strongest_dependency()[2])
                except AssertionError:
                    raise
                except KeyboardInterrupt:
                    raise
                except NameError:
                    raise
                except:
                    print 'Problem testing %s' % matrix
                    #print sys.exc_info()[0]
                    print sys.exc_info()[1]
                    #raise
            print 'Mean LL ratio: ', N.mean(LL_ratios)
            pylab.figure()
            pylab.scatter(best_mis, LL_ratios)
            pylab.figure()
            pylab.scatter(best_mis, num_seqs)
            #raise ''

    # test using custom data
    if True:
        num_features = 4
        num_values = 2
        data = [
          [0,0,1,1],
          [0,0,1,0],
          [1,0,1,1],
          [1,0,1,0],
          [0,0,1,1],
          [0,0,1,0],
          [1,0,1,1],
          [1,0,1,0],
          [0,1,0,0],
          [0,1,0,1],
          [1,1,0,0],
          [1,1,0,1],
          [0,1,0,0],
          [0,1,0,1],
          [1,1,0,0],
          [1,1,0,1],
        ]
        print cross_validate(data, num_features, num_values) / len(data)
        dependencies = DependencyAnalyser(data, num_features, num_values)
        dependencies.highlight_tree_edges()
        dependencies.write_graph('graph')

        #for datum in data:
        #  print datum, calculate_log_likelihood_ratio(datum, conditioners, log_marginals, log_conditionals)
