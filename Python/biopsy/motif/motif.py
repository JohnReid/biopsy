#
# Copyright John Reid 2007
#

import numpy, math
import biopsy

#mutual_informations = []

class Motif( object ):
    """
    Takes a sequence of sequences and builds joint and marginal probabilities
    Also calculates mutual information between all pairs of bases
    """

    def __init__( self, seqs, marginal_prior = 0.0, joint_prior = 0.0, mi_threshold = 0.0 ):
        self.seqs = seqs
        self.N = len(seqs)
        self.K = self.N and len(seqs[0]) or 0
        self.calculate_joint_and_marginals(marginal_prior, joint_prior)
        self.build_tree(mi_threshold)

    def normalise_joint_and_marginals(self):
        for i in xrange(self.K):
            self.marginal[i] = normalise(self.marginal[i])
            for j in xrange(self.K):
                self.joint[i,j] = normalise(self.joint[i,j])

    def calculate_joint_and_marginals(self, marginal_prior, joint_prior):
        self.marginal = numpy.zeros( (self.K,4) ) + marginal_prior
        self.joint = numpy.zeros( (self.K,self.K,4,4) ) + joint_prior
        for s in self.seqs:
            assert self.K == len(s), 'All sequences should be of same length'
            for i,bi in enumerate(s):
                if bi != 4:
                    self.marginal[i,bi] += 1.0
                    for j,bj in enumerate(s):
                        if bj != 4:
                            self.joint[i,j,bi,bj] += 1.0
        self.normalise_joint_and_marginals()
        self.mutual_information = numpy.zeros( (self.K,self.K) )
        for i in xrange(self.K):
            for bi in xrange(4):
                marginal_i = self.marginal[i,bi]
                if marginal_i > 0.0:
                    for j in xrange(self.K):
                        for bj in xrange(4):
                            marginal_j = self.marginal[j,bj]
                            if marginal_j > 0.0:
                                joint = self.joint[i,j,bi,bj]
                                if joint > 0.0:
                                    self.mutual_information[i,j] += joint * (
                                            math.log(joint)
                                            - math.log(marginal_i)
                                            - math.log(marginal_j)
                                    )
                                    #mutual_informations.append(self.mutual_information[i,j])

    def build_tree(self, mi_threshold = 0.0):
        """
        Builds a tree from the mutual information between the base distributions
        """
        import boost.graph as bgl
        # first build an (undirected) graph of the (negative) mutual informations
        g = bgl.Graph()
        vertices = [ g.add_vertex() for i in xrange(self.K) ]
        weight_map = g.add_edge_property('weight', 'float')
        joint_probabilities = g.add_edge_property('joints')
        marginal_probabilities = g.add_vertex_property('marginals')

        # add the edges and marginal probs
        for i in xrange(self.K):
            marginal_probabilities[vertices[i]] = self.marginal[i]

        # add the edges and conditional probs
        for i in xrange(self.K):
            for j in xrange(i+1,self.K):
                if self.mutual_information[i,j] > mi_threshold:
                    e = g.add_edge(vertices[i],vertices[j])
                    weight_map[e] = -self.mutual_information[i,j]
                    joint_probabilities[e] = self.joint[i,j]

        # find the minimum spanning tree and remove those edges not in it
        spanning_tree = bgl.kruskal_minimum_spanning_tree(g,weight_map)
        for e in g.edges:
            if e not in spanning_tree:
                g.remove_edge(e)

        def index_of(v):
            for i in xrange(self.K):
                if vertices[i] == v:
                    return i
            raise RuntimeError( 'Could not find vertex: ' + v )

        # build a directed forest from the (undirected) graph
        def build_forest(g, root_idx = None):
            "Take a forest and build an isomorphic Digraph"
            import boost.graph as bgl

            vertices_visited = set()
            diforest = bgl.Digraph()
            edge_label = diforest.add_edge_property( 'label', 'string' )
            di_conditionals = diforest.add_edge_property( 'conditionals' )
            di_marginals = diforest.add_vertex_property( 'marginals' )
            di_vertices = [ diforest.add_vertex() for i in xrange(self.K) ]

            def visit(i):
                if i in vertices_visited: return
                vertices_visited.add(i)
                v = vertices[i]
                di_v = di_vertices[i]
                di_marginals[di_v] = marginal_probabilities[v]
                for neighbour in g.adjacent_vertices(v):
                    j = index_of(neighbour)
                    if j not in vertices_visited:
                        e = diforest.add_edge(di_v, di_vertices[j])
                        edge_label[e] = '%.3f' % (-weight_map[g.edge(v,vertices[j])])
                        joint = joint_probabilities[g.edge(v,neighbour)]
                        if index_of(neighbour) < index_of(v): joint = joint.T # make sure we have correct orientation of joint
                        di_conditionals[e] = joint_to_conditional(joint)
                        visit(j) # recurse

            # visit all the vertices - starting with root if given
            if None != root_idx: visit(root_idx)
            for i in xrange(self.K): visit(i)

            # check each node has max one parent
            for v in diforest.vertices: assert diforest.in_degree(v) < 2

            return diforest

        # set the root to have the highest MI
        root_idx = self.mutual_information.argmax() % self.mutual_information.shape[0]
        self.d = build_forest(g, root_idx = root_idx)

    def log_likelihood_from_pssm(self, seq):
        seq = as_seq(seq)
        assert len(seq) == len(self.marginal), 'Sequence and pssm must be same length'
        log_pssm = numpy.log(self.marginal)
        return sum( [ l[b] for l,b in zip(log_pssm,seq) ] )

    def log_likelihood_from_forest(self, seq):
        vertices = [ v for v in self.d.vertices ]
        assert len(seq) == len(vertices), 'Sequence and forest must be same length: %d, %d' % (len(seq), len(vertices))
        marginals = self.d.vertex_properties['marginals']
        conditionals = self.d.edge_properties['conditionals']
        def index_of(v): return vertices.index(v)
        def LL_for(v,b):
            if 0 == self.d.in_degree(v):
                p = marginals[v][b]
                return numpy.log(p)
            elif 1 == self.d.in_degree(v):
                e = self.d.in_edges(v).next()
                u = self.d.source(e)
                i = index_of(u)
                assert b == seq[index_of(v)]
                p = conditionals[e][seq[i],b]
                return numpy.log(p)
            else: raise RuntimeError('Can only handle vertices with at most one parent')
        return sum( [ LL_for(v,b) for v,b in zip(vertices, seq) ] )

    def LL_for_seqs(self, seqs, method):
        return sum( [ self.method(s) for s in seqs ] )

    def pssm_LL_for_seqs(self, seqs): return self.LL_for_seqs(seqs, Motif.log_likelihood_from_pssm)
    def forest_LL_for_seqs(self, seqs): return self.LL_for_seqs(seqs, Motif.log_likelihood_from_forest)

    def LL_ratio_for_seqs(self, seqs):
        return sum( [ self.log_likelihood_from_forest(s) - self.log_likelihood_from_pssm(s) for s in seqs ] )

def normalise(a):
    s = a.sum()
    if 0.0 == s: return a
    else: return a / s

def joint_to_conditional(A):
    "Returns a conditional distribtion from a joint one"
    return numpy.array( [ normalise(a) for a in A ] )

def index_from_base(b, allow_unknowns = True):
    if 'a' == b or 'A' == b: return 0
    elif 'c' == b or 'C' == b: return 1
    elif 'g' == b or 'G' == b: return 2
    elif 't' == b or 'T' == b: return 3
    elif allow_unknowns and ('n' == b or 'N' == b): return 4
    else: raise RuntimeError( 'Unknown base: ' + str(b) )

def seq_from_string(s, allow_unknowns = True):
    return numpy.array( [ index_from_base(b, allow_unknowns) for b in s ], dtype=numpy.int32 )

def as_seq(s):
    if isinstance(s, str): return(seq_from_string(s))
    else: return numpy.asarray(s, dtype = numpy.int32)

def seqs_from_strings(strings, allow_unknowns = True):
    return numpy.array( [ seq_from_string(s, allow_unknowns) for s in strings ] )

def write_dependencies(g, basename):
    "Take the dependency graph g and layout dependencies as SVG graph"
    dot_filename = '%s.dot' % basename
    svg_filename = '%s.svg' % basename

    from copy import copy
    g_to_write = copy(g)
    del g_to_write.edge_properties['conditionals']
    del g_to_write.edge_properties['object']
    del g_to_write.vertex_properties['marginals']
    del g_to_write.vertex_properties['object']

    g_to_write.write_graphviz(dot_filename)
    import os
    os.system('neato -Elen=2 -Tsvg %s -o%s' % (dot_filename, svg_filename))


def transfac_matrix_sequences():
    import biopsy.transfac as t
    for m in t.Matrix.all():
        seqs = m.sequences
        if len(seqs):
            yield m, seqs


class Parameters(object):
    def __init__(self, marginal_prior = 0.0, joint_prior = 0.0, mi_threshold = 0.0):
        self.marginal_prior = marginal_prior
        self.joint_prior = joint_prior
        self.mi_threshold = mi_threshold

    def build_motif(self, seqs):
        return Motif(seqs, self.marginal_prior, self.joint_prior, self.mi_threshold)

    def k_fold_cross_validation_test(self, seqs, K = 5):
        from infpy import k_fold_cross_validation
        LL_ratio = 0.0
        for training, validation in k_fold_cross_validation(seqs, K, randomise = True):
            if not len(validation): continue
            motif = self.build_motif(training)
            LL_ratio += motif.LL_ratio_for_seqs(validation)/len(validation)
        return LL_ratio/(K*len(motif.marginal))

    def __str__(self):
        return 'marginal_prior=%.2f, joint_prior=%.4f, mi_threshold=%.2f' %(
    self.marginal_prior,
    self.joint_prior,
    self.mi_threshold
        )

def sequences_from_jaspar_file(jaspar_file):
    length = None
    for l in jaspar_file:
        if l.startswith('>'): continue
        if not l.strip(): continue
        site = ''.join( [ c for c in l if c.isupper() ] )
        if length == None: length = len(site)
        elif len(site) != length: continue
        yield site

_jaspar_dir = os.path.join(biopsy.get_data_dir(), 'Jaspar', 'JASPAR_CORE')
_jaspar_phylofacts_dir = os.path.join(biopsy.get_data_dir(), 'Jaspar', 'JASPAR_PHYLOFACTS')
def jaspar_sequences(dir = _jaspar_dir):
    import os
    for filename in os.listdir(dir):
        f = open(os.path.join(dir, filename), 'r')
        yield (
                filename.split('.')[0],
                [ s for s in sequences_from_jaspar_file(f) ]
        )
def jaspar_phylofacts_sequences():
    for x in jaspar_sequences(_jaspar_phylofacts_dir):
        yield x

if '__main__' == __name__:

    class LL_Ratio_Tests(object):
        parameters = [
                Parameters(.25,0.0625,0.0),
                Parameters(.25,0.0625,0.2),
                Parameters(.25,0.0625,0.4),
                Parameters(.25,0.0625,0.6),
                Parameters(.25,0.0625,0.8),
                #Parameters(1,1,0.0),
                #Parameters(1,1,0.2),
                #Parameters(1,1,0.4),
                #Parameters(1,1,0.6),
                #Parameters(1,1,0.8),
        ]

        from itertools import islice
        sequences = {
                'transfac' : (list(transfac_matrix_sequences()),'b'),
                'jaspar core' : (list(jaspar_sequences()),'r'),
                #'jaspar phylofacts' : (islice(jaspar_phylofacts_sequences(),4),'y'),
        }

        def test_params(self):
            from itertools import islice
            import pylab
            for i, params in enumerate(self.parameters):
                test_name = str(params)
                pylab.figure()
                pylab.title(test_name)
                for seqs_name, (sequences, c) in self.sequences.iteritems():
                    lengths = []
                    ratios = []
                    LL_ratio_total = 0.0
                    num_sources = 0
                    for m, seqs in sequences:
                        ratio = params.k_fold_cross_validation_test(seqs_from_strings(seqs))
                        LL_ratio_total += ratio
                        lengths.append(len(seqs))
                        ratios.append(ratio)
                        num_sources += 1
                    print 'average LL_ratio: %.3f : %s' % (LL_ratio_total/num_sources, test_name)
                    pylab.scatter(lengths, ratios, label=seqs_name, c=c)
                    pylab.xlabel('# seqs')
                    pylab.ylabel('LL ratio/per base')
                pylab.xlim(xmin=0)
                pylab.axhspan(0,0) # draw a line through LL_ratio = 0
                pylab.legend
                pylab.savefig('motif_%d.png'%i)
                pylab.savefig('motif_%d.ps'%i)

    LL_Ratio_Tests().test_params()

    if False:
        class MotifTester(object):
            def __init__(self, seqs, name = 'tree'):
                self.seqs = seqs_from_strings(seqs)
                self.motif = Motif(self.seqs)
                self.g = self.motif.build_tree()
                write_dependencies(self.g, 'tree')

            def __call__(self, test_seq):
                test_seq
                print 'LL from pssm:   %.5f : %s' % (self.motif.log_likelihood_from_pssm(seq_from_string(test_seq)), test_seq)
                print 'LL from forest: %.5f : %s' % (self.motif.log_likelihood_from_forest(seq_from_string(test_seq)), test_seq)
                print

        t = MotifTester(
                [
                  'acgtacgt',
                  'tcgtacgt',
                  'acctcagt',
                  'agctcagt',
                ]
        )
        t('agctcagt')
        t('tgctcagt')

        t = MotifTester(
                [
                  'at',
                  'at',
                  'cg',
                  'cg',
                ]
        )
        t('at')
        t('ag')
