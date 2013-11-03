#
# Copyright John Reid 2007
#

print 'testing gapped pssms'

import numpy

try: hmm
except: import biopsy.gapped_pssms as hmm


def print_model( m ):
    print 'rho'
    print m.var_dist.rho
    print 'eta'
    print m.var_dist.eta
    print 'tau'
    print m.var_dist.tau


def create_model(
        K,
        sequences,
        psi = 0.01 * numpy.ones( 4 ),
        theta = 100.0 * numpy.ones( 4 ),
        phi = numpy.array( [ 9.5, 0.5 ] ),
        upsilon = numpy.array( [ 99.0, 1.0 ] )
):
    print 'testing state map'
    state_map = hmm.StateMap( K )
    hmm.test_state_map( state_map )
    assert K == state_map.K
    state_map.S
    state_map.b
    state_map.c
    state_map.g
    state_map.k
    state_map.m
    state_map.s

    print 'Creating sequences'
    s = hmm.ObservedSequences( sequences )
    print s.sequences
    assert s.N == len( sequences )
    assert s.I(0) == len( sequences[ 0 ] )

    print 'Creating data'
    d = hmm.ObservedData(
            K,
            s,
            psi,
            theta,
            phi,
            upsilon
    )
    assert psi.all() == d.Psi.all()
    assert theta.all() == d.Theta.all()
    assert phi.all() == d.Phi.all()
    assert upsilon.all() == d.Upsilon.all()

    print 'Drawing hidden data'
    hd = hmm.HiddenData( d )
    print hd.draw_sequence( 10 )

    print 'Creating model'
    m = hmm.Model( d )
    m.data
    m.predecessor_states
    m.var_dist
    # print_model( m )

    return m


psi = 0.01 * numpy.ones( 4 )
theta = 100.0 * numpy.ones( 4 )
phi = numpy.array( [ 0.01, 0.01 ] )
upsilon = numpy.array( [ 0.01, 0.01 ] )


def test_transitions():
    "Test that transition probs are set correctly. Can check by inspecting svg."
    K = 3
    sequences =     [
            'tgacg',
    ]
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    v = m.var_dist
    v.tau = numpy.array(
            [
                    [ .90, .10 ],
                    [ .50, .50 ],
                    [ .80, .20 ]
            ]
    )
    hmm.write_model_svg(
            m,
            name = 'transitions',
            dir = 'test',
            show_rev_comp = True,
            show_dists = False,
            edge_lengths = 3.0
    )
    return m, v
# m, v = test_transitions()

def test_eta_simple():
    "Test update eta with simple 1 base pssm"
    K = 1
    sequences =     [
            'tg',
    ]
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    v = m.var_dist
    v.rho = numpy.array(
            [
                    [
                            [ .98, .01, .01 ],
                            [ .01, .98, .01 ]
                    ]
            ]
    )
    print 'Updating'
    for i in xrange( 100 ):
        m.update(
                update_rho = False,
                update_eta = True,
                update_tau = False
        )
    hmm.write_model_svg(
            m,
            name = 'eta_simple',
            dir = 'test',
            show_rev_comp = True,
            show_dists = True,
            edge_lengths = 3.0
    )
    return m, v
# m, v = test_eta_simple()

def test_eta_gapped():
    "Test update eta with gapped 2 base pssm"
    K = 2
    sequences =     [
            'tag', # 'a' is the gap
            #'tg',
    ]
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    v = m.var_dist
    v.rho = [
            numpy.array(
                    [
                            [ .01, .94, .01, .01, .01, .01, .01 ],
                            [ .01, .01, .94, .01, .01, .01, .01 ],
                            [ .01, .01, .01, .94, .01, .01, .01 ],
                    ]
            ),
#               numpy.array(
#                       [
#                               [ .01, .94, .01, .01, .01, .01, .01 ],
#                               [ .01, .01, .01, .94, .01, .01, .01 ],
#                       ]
#               )
    ]
    for i in xrange( 100 ):
        m.update(
                update_rho = False,
                update_eta = True,
                update_tau = False
        )
    print 'Updating'
    hmm.write_model_svg(
            m,
            name = 'eta_gapped',
            dir = 'test',
            show_rev_comp = True,
            show_dists = True,
            edge_lengths = 3.0
    )
    return m, v
# m, v = test_eta_gapped()

def test_tau_simple():
    K = 1
    sequences =     [
            'tg',
    ]
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    v = m.var_dist
    v.rho = numpy.array(
            [
                    [
                            [ .98, .01, .01 ],
                            [ .01, .98, .01 ]
                    ]
            ]
    )
    print 'Updating'
    for i in xrange( 100 ):
        m.update(
                update_rho = False,
                update_eta = True,
                update_tau = True
        )
    hmm.write_model_svg(
            m,
            name = 'tau_simple',
            dir = 'test',
            show_rev_comp = True,
            show_dists = False,
            edge_lengths = 3.0
    )
    return m, v
# m, v = test_tau_simple()




def test_tau_gapped():
    K = 2
    sequences =     [
            'tag',
            'tg',
    ]
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    v = m.var_dist
    v.rho = [
            numpy.array(
                    [
                            [ .01, .94, .01, .01, .01, .01, .01 ],
                            [ .01, .01, .94, .01, .01, .01, .01 ],
                            [ .01, .01, .01, .94, .01, .01, .01 ],
                    ]
            ),
            numpy.array(
                    [
                            [ .01, .94, .01, .01, .01, .01, .01 ],
                            [ .01, .01, .01, .94, .01, .01, .01 ],
                    ]
            )
    ]
    print 'Updating'
    for i in xrange( 100 ):
        m.update(
                update_rho = False,
                update_eta = True,
                update_tau = True
        )
    hmm.write_model_svg(
            m,
            name = 'tau_gapped',
            dir = 'test',
            show_rev_comp = True,
            show_dists = False,
            edge_lengths = 3.0
    )
    return m, v
# m, v = test_tau_gapped()


def test_rho_simple():
    K = 2
    sequences = [
            'tga',
    ]
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    v = m.var_dist
    v.eta = numpy.array( # a pssm that looks like G[G]A
            [
                    [ 0.25, 0.25, 0.25, 0.25 ],
                    [ 0.01, 0.01, 0.97, 0.01 ],
                    [ 0.01, 0.01, 0.97, 0.01 ],
                    [ 0.97, 0.01, 0.01, 0.01 ]
            ]
    )
    print 'Updating'
    for i in xrange( 100 ):
        print m.update(
                update_rho = True,
                update_eta = False,
                update_tau = False
        ),
        print "\n".join( str( r ) for r in v.r_mode )
    hmm.write_model_svg(
            m,
            name = 'rho_simple',
            dir = 'test',
            show_rev_comp = True,
            show_dists = False,
            edge_lengths = 3.0
    )
    return m, v
m, v = test_rho_simple()




def test_simple():
    K = 2
    sequences = [
            'tg',
            'tg',
            'tag',
            'tag',
    ]
    print theta
    m = create_model( K, sequences, psi, theta, phi, upsilon )
    #hmm.write_model_svg(
    #       m,
    #       name = '%d' % K,
    #       dir = 'test',
    #       show_rev_comp = True,
    #       show_dists = False,
    #       edge_lengths = 3.0
    #)
    #raise ""
    v = m.var_dist
    v.eta = numpy.array( # a pssm that looks like T[A]G
            [
                    [ 0.25, 0.25, 0.25, 0.25 ],
                    [ 0.01, 0.01, 0.01, 0.97 ],
                    [ 0.97, 0.01, 0.01, 0.01 ],
                    [ 0.01, 0.01, 0.97, 0.01 ]
            ]
    )
    print 'Updating'
    for i in xrange( 100 ):
        m.update(
                update_rho = True,
                update_eta = False,
                update_tau = True
        )
    print "\n".join(str( r ) for r in m.var_dist.r_mode)
    hmm.write_model_svg(
            m,
            name = 'simple',
            dir = 'test',
            show_rev_comp = True,
            show_dists = True,
            edge_lengths = 3.0
    )
    return m, v
# m, v = test_simple()






def test_write_models():
    for k in xrange( 10 ):
        model = hmm.Model(
                hmm.ObservedData(
                        k + 1,
                        hmm.ObservedSequences( sequences )
                )
        )
        hmm.write_model_svg(
                model,
                'states_%d' % (k + 1),
                show_rev_comp = False
        )
# test_write_models()
