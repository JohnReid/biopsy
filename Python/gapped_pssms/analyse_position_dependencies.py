#
# Copyright John Reid 2008,2009
#

import sys, os, logging, itertools
from hmm.pssm import seq_to_numpy, numpy_to_seq
import hmm.pssm.logo as L
from biopsy.chow_liu.chow_liu import DependencyAnalyser
from gapped_pssms.parse_gapped_format import parse_models, build_hmm_from_semi_parsed
from cookbook.dicts import DictOf
from cookbook.pylab_utils import pylab_ioff
import pylab as P
import numpy


@pylab_ioff
def conditional_logo(joint):
    """
    Produces an image summarising the marginal, conditional and joint distributions of 2 bases.
    """
    import pylab as P, numpy as N, hmm.pssm.logo as L

    #dpi = 150
    #fig = P.figure(figsize=(9,10), dpi=dpi, facecolor='white')

    #dpi = 150
    fig = P.figure(figsize=(6, 6), facecolor='white')
    dpi = fig.dpi

    def logo_from_dist(dist):
        return L.dist_as_image(dist/dist.sum(), (dpi, dpi))

    def place_logo(logo, x, y):
        P.figimage(N.asarray(logo, dtype=N.float32) / 255., x*dpi, y*dpi)

    # conditionals: logo
    for x in xrange(4):
        cond = joint[x]
        place_logo(logo_from_dist(cond/cond.sum()), x+1, 0)
    for y in xrange(4):
        cond = joint[:,y]
        place_logo(logo_from_dist(cond/cond.sum()), 5, 4-y)

    # marginals: logo
    x_marg = joint.sum(axis=1)
    place_logo(logo_from_dist(x_marg/x_marg.sum()), 2.5, 5)
    y_marg = joint.sum(axis=0)
    place_logo(logo_from_dist(y_marg/y_marg.sum()), 0, 2.5)

    # joint distribution: heat map
    Z = N.ones((dpi, dpi))
    for x in xrange(4):
        for y in xrange(4):
            area = joint[x,y]
            edge_len = int(area * dpi)
            if edge_len:
                offset = (dpi-edge_len)/2.
                P.figimage(N.ones((edge_len, edge_len)), xo=dpi*(x+1)+offset, yo=dpi*(4-y)+offset)

    return fig
#logo = conditional_logo(joint)
#logo.savefig('dependencies/test.png', dpi=logo.dpi)
#raise ''



@pylab_ioff
def conditional_logo_2(joint, X=None, Y=None):
    """
    Produces an image summarising the marginal, conditional and joint distributions of 2 bases.
    """
    import pylab as P, numpy as N, hmm.pssm.logo as L
    from PIL import Image, ImageDraw

    #dpi = 150
    #fig = P.figure(figsize=(9,10), dpi=dpi, facecolor='white')

    #dpi = 150
    if None == X:
        X = joint.shape[0]
    if None == Y:
        Y = joint.shape[1]
    fig = P.figure(figsize=(X+2, Y+2), facecolor='white')
    dpi = fig.dpi

    def place_image(logo, x, y):
        P.figimage(N.asarray(logo, dtype=N.float32) / 255., x*dpi, y*dpi)

    def place_square(area, x, y, colour=N.array([255,255,255])):
        edge_len = int((area**.5)*dpi)
        if edge_len:
            offset = (dpi-edge_len)/2.
            array = N.ones((edge_len, edge_len, 3))
            array[:,:] = colour
            P.figimage(array, xo=dpi*x+offset, yo=dpi*y+offset)

    # base labels
    bases = ['A', 'C', 'G', 'T', None]
    font = L.get_font(L._default_font, font_size=int(dpi/2))
    def image_for_base(b):
        image = Image.new('RGB', (dpi, dpi), 'white')
        draw = ImageDraw.Draw(image)
        if b:
            textsize = draw.textsize(b, font=font)
            draw.text((dpi/2-textsize[0]/2,dpi/2-textsize[1]/2), b, font=font, fill=colour)
        else:
            draw.rectangle(((dpi/4,dpi/4), (3*dpi/4,3*dpi/4)), outline=colour)
        return image

    #colour = '#009999'
    colour = 'black'
    for x, image in enumerate(map(image_for_base, bases[:X])):
        place_image(image, x+1, Y+1)

    #colour = '#990099'
    colour = 'black'
    for y, image in enumerate(map(image_for_base, bases[:Y])):
        place_image(image, 0, Y-y)

    # joint distribution
    for x in xrange(X):
        for y in xrange(Y):
            place_square(joint[x,y], x+1, Y-y, colour=N.array([128,128,128]))

    # marginals
    x_marg = joint.sum(axis=1)
    for x in xrange(X):
        place_square(x_marg[x], x+1, 0)
    y_marg = joint.sum(axis=0)
    for y in xrange(Y):
        place_square(y_marg[y], X+1, Y-y)

    return fig
#logo = conditional_logo_2(joint, X=5, Y=4)
#logo.savefig('dependencies/test.png', dpi=logo.dpi)
#raise ''




def sequences_from_fasta(fasta):
    """Yields sequences from fasta file."""
    import corebio.seq_io.fasta_io
    return itertools.imap(
      lambda s: s.strip('nN'),
      itertools.imap(
        str,
        corebio.seq_io.fasta_io.iterseq(
          open(fasta, 'r'),
          corebio.seq.dna_alphabet
        )
      )
    )


def sites_from_states(states, sequence, background_states):
    assert len(states) == len(sequence)
    site_seq, site_states = [], []
    in_site = False
    for state, base in zip(states, sequence):
        if in_site and state in background_states:
            yield site_seq, site_states
            site_seq, site_states = [], []
        in_site = state not in background_states
        if in_site:
            site_seq.append(base)
            site_states.append(state)
    if in_site:
        yield site_seq, site_states


def uncomplement_base(base, state, traits):
    """
    Takes the base and its state and if it is a reverse complement state returns its complement.
    """
    if state in traits.reverse_complements:
        return 3-base, traits.reverse_complements[state]
    else:
        return base, state

def uncomplement_site(site_seq, site_states, traits):
    """
    Un-complements a site.
    """
    seq, states = [], []
    for base, state in zip(site_seq, site_states):
        base, state = uncomplement_base(base, state, traits)
        seq.append(base)
        states.append(state)
    return seq, states



def fisher_test(table, B=2000):
    import rpy2.robjects.numpy2ri
    from rpy2.robjects import r
    args = {
        'simulate.p.value' : True,
        'B'                : B
    }
    return r['fisher.test'](table, **args)


#
# Initialise the logging
#
logging.basicConfig(level=logging.INFO)
file_handler = logging.FileHandler('position-dependencies.log')
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger('').addHandler(file_handler)


P.rcParams['figure.dpi'] = 150
bases = ['a', 'c', 'g', 't']
p_binding_site = 0.001
min_mutual_info = 0.3
ignore_edges_below = .1
image_size = (200,500)
model_dir = sys.argv[1]
sequence_dir = sys.argv[2]
dependency_dir = os.path.join(model_dir, 'dependencies')
sequence_filename_fmt = '%strimRM.fa'
min_num_sites = 20

def pssms():
    for pssm in [
        'T99002',
        'T99003',
        'T99004',
        'T99005',
        'T99006',
    ]:
        for i in xrange(10):
            yield pssm, '%03d' % i

#def pssms():
#    yield 'T99006', '009'

#pssms = [
#  ('T99002', '000'),
#  ('T99003', '000'),
#  ('T99004', '000'),
#  ('T99005', '000'),
#  ('T99006', '000'),
#]
#model_dir = os.path.join('c:\\', 'Johns', 'Writing', 'GappedPssms', 'Single-Gap', 'results-2')
#model_dir = '/home/reid/Analysis/GappedPssms/apr-2009/single-gap'
#sequence_dir = '/home/reid/Data/GappedPssms/apr-2009/'


fisher_p_values = list()
for fragment, pssm in pssms():
    sequence_file = os.path.join(sequence_dir, sequence_filename_fmt % fragment)
    model_file = os.path.join(model_dir, '%s-%s.pssm' % (fragment, pssm))

    logging.info('Loading sequences: %s', sequence_file)
    sequences = list(sequences_from_fasta(sequence_file))
    numpy_seqs = map(seq_to_numpy, sequences)
    logging.info('Loaded %d sequences', len(sequences))


    logging.info('Parsing PSSMs: %s', model_file)
    pssms = list(parse_models(open(model_file)))


    logging.info('Building models')
    models = [
      build_hmm_from_semi_parsed(parsed, p_binding_site=p_binding_site)
      for parsed in pssms
    ]

    def nucleotide_dist():
        return numpy.zeros(4) + .25
    base_dists = DictOf(nucleotide_dist)

    min_site_length = 20
    logging.info('Analysing sequences')
    for hmm, traits in models:
        sites = []
        for sequence in numpy_seqs:

            # analyse the sequence for its most likely state sequence
            LL, states = hmm.viterbi(sequence)

            # for each site
            for site_seq, site_states in sites_from_states(states, sequence, traits.background_states):

        # is it long enough?
                if len(site_seq) > min_site_length:
                    # store uncomplemented version of site
                    sites.append(uncomplement_site(site_seq, site_states, traits))
        logging.info('Found %d sites', len(sites))
        if len(sites) < min_num_sites:
            logging.info('Not enough sites')
            continue

        # work out which states are in the sites
        states_in_sites = set()
        gap_state = None
        for site_seq, site_states in sites:
            states_in_sites.update(site_states)
        # find the gap base, it is the one with an even state index
        for state in states_in_sites:
            if 0 == state % 2:
                assert None == gap_state
                gap_state = state
        states_in_sites = list(states_in_sites)
        states_in_sites.sort()
        state_to_feature = dict((s, f) for f, s in enumerate(states_in_sites))
        num_features = len(state_to_feature)
        logging.info('Gap state=%d; feature index=%d', gap_state, state_to_feature[gap_state])

        # transform sites into format suitable for Chow-Liu analysis
        sites_as_features = []
        for site_seq, site_states in sites:
            site_features = [4] * num_features
            for base, state in zip(site_seq, site_states):
                site_features[state_to_feature[state]] = base
                if 4 != base:
                    base_dists[state_to_feature[state]][base] += 1
            sites_as_features.append(site_features)

                #print site_seq
                #print site_states
            #print 'Done sequence'

        dependencies = DependencyAnalyser(
            sites_as_features,
            num_features,
            5,
            pseudo_count=0.0,
            min_mutual_info=min_mutual_info
        )
        dependencies.remove_edges_below(ignore_edges_below)
        dependencies.highlight_tree_edges()
        dependency_graph_name = os.path.join(dependency_dir, 'dependencies-%s-%s' % (fragment, pssm))
        logging.info('Saving dependency graph to %s', dependency_graph_name)
        dependencies.write_graph(dependency_graph_name)
        gap_feature = state_to_feature[gap_state]
        for (i1, i2), counts in dependencies.counts.iteritems():
            if gap_feature in (i1, i2):
                if i1 == i2:
                    continue
                elif gap_feature == i1:
                    counts = counts[:,:4]
                elif gap_feature == i2:
                    counts = counts[:4,:]
                counts = counts.take(numpy.where(counts.sum(axis=1) > 0)[0], axis=0)
                counts = counts.take(numpy.where(counts.sum(axis=0) > 0)[0], axis=1)
                if counts.shape[0] > 1 and counts.shape[1] > 1:
                    fisher_test_results = fisher_test(counts)
                    fisher_p_value = fisher_test_results[0][0]
                    fisher_p_values.append(fisher_p_value)
                    if fisher_p_value < 1e-4:
                        logging.info('%s-%s: Fisher test: %2d %2d %e', fragment, pssm, i1, i2, fisher_p_value)

        most_dependent_pair = dependencies.strongest_dependency()[0]
        #most_dependent_pair = None
        gap_feature = state_to_feature[gap_state]
        for pair, joint, mi in dependencies.mutual_infos:
            if mi > min_mutual_info and (pair == most_dependent_pair or gap_feature in pair):
                logging.info('%d-%d: mutual info=%.3f', pair[0], pair[1], mi)
                cond_logo = conditional_logo_2(joint, X=4+int(pair[0]==gap_feature), Y=4+int(pair[1]==gap_feature))
                image_filename = os.path.join(dependency_dir, 'dependencies-%s-%s-%02dx%02d.png' % (fragment, pssm, pair[1], pair[0]))
                cond_logo.savefig(image_filename, dpi=cond_logo.dpi)
