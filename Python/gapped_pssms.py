#
# Copyright John Reid 2006
#

from biopsy.gapped_pssms import *

def multi_model_learn(
        seqs,
        K,
        alpha,
        phi,
        varphi,
        num_models = 5
):
    models = [
            VariationalModel(
                    seqs = seqs,
                    K = K,
                    alpha = alpha,
                    phi = phi,
                    varphi = varphi
            )
            for i in xrange( num_models )
    ]
    learners = [ model.learner() for model in models ]
    for i in xrange( 12 ):
        for i, learner in enumerate( learners ):
            if None == learner.gain or learner.gain > 1e-5:
                learner()
        yield learners

def test_multi_model(
        seqs,
        K,
        phi,
        varphi,
        alpha
):
    import pylab
    for learners in multi_model_learn(
            seqs = seqs,
            K = K,
            alpha = alpha,
            phi = phi,
            varphi = varphi,
            num_models = 10
    ):
        pylab.clf()
        for learner in learners:
            pylab.plot( learner.LL )

def test_entropy(
        seqs,
        K,
        phi,
        varphi,
        alpha
):
    var_model = VariationalModel(
            seqs = seqs,
            K = K,
            alpha = alpha,
            phi = phi,
            varphi = varphi
    )
    import pylab
    pssm_entropy = []
    start_entropy = []
    has_gap_entropy = []
    bg_gap_entropy = []
    for i in xrange( 6 ):
        var_model.update()
        pssm_entropy.append( var_model.pssm_entropy_per_base() )
        start_entropy.append( var_model.start_entropy_per_sequence() )
        has_gap_entropy.append( var_model.has_gap_entropy_per_sequence() )
        bg_gap_entropy.append( var_model.bg_entropy() )
        pylab.clf()
        pylab.plot( pssm_entropy, label='pssm entropy' )
        pylab.plot( start_entropy, label='start entropy' )
        pylab.plot( has_gap_entropy, label='has gap entropy' )
        pylab.plot( bg_gap_entropy, label='background entropy' )
        pylab.legend()

def base_to_index( b ):
    "Convert a,c,g,t to 0,1,2,3"
    if 'a' == b or 'A' == b: return 0
    elif 'c' == b or 'C' == b: return 1
    elif 'g' == b or 'G' == b: return 2
    elif 't' == b or 'T' == b: return 3
    else:
        raise RuntimeError( 'Unknown base "%s"' % b )

def index_to_base( i ):
    if 0 == i: return 'a'
    if 1 == i: return 'c'
    if 2 == i: return 'g'
    if 3 == i: return 't'
    if 4 == i: return 'n'
    else:
        raise RuntimeError( 'Unknown index: %d' % i )

def seq_to_numeric( seq ):
    "Convert a sequence of a,c,g,t to a sequence of 0,1,2,3"
    return [ base_to_index( b ) for b in seq ]

def create_model( K, seqs, c_impl = True ):
    if c_impl:
        model = VariationalModel_C.create( K, seqs )
    else:
        model = VariationalModel(
                seqs = seqs,
                K = K,
                alpha = [ 1.0 ] * 2,
                phi = [ 0.1 ] * K,
                varphi = [ 10.0 ] * K
        )
    return model

def learn_model(
        model,
        logo_file = None,
        max_updates = 25
):
    learner = VariationalModel.Learner(
            model,
            max_updates = max_updates
    )
    learner.learn( verbose = True )
    if None != logo_file:
        format_weblogo_from_dist(
                model.omega_array()[1:,],
                logo_file,
                'png',
        )

def simple_test():
    """A simple test for gapped pssms.

    We use a motif AAAA with a gap at position 2.
    """
    seqs = [
            'gaaaac',
            'aacaat',
            'aacaat',
            'aacaag',
            'aaaatc',
            'anaacg',
    ]
    model = create_model( 4, seqs )
    learn_model( model, 'pssms/simple_test' )
    return model


def discrete_KL( p, q ):
    return sum(
            [
                    _p * ( math.log(_p) - math.log(_q) )
                    for _p, _q
                    in zip( p, q )
            ]
    )

def compare_model_with_truth( model, truth ):

    print 'Comparing learnt model with truth....'

    locations_difference = [
            l - numpy.argmax( model.nu_sequence()[i] )
            for i, l
            in enumerate( truth.locations )
    ]
    print 'Site location offsets:', locations_difference

    assert len( truth.has_gap ) == len( model.mu_array() )
    gap_correct = [
            t == ( g > 0.5 )
            for t, g
            in zip( truth.has_gap, model.mu_array() )
    ]
    print 'Correct gap predictons:', gap_correct

    print 'Gap position offset:', truth.pssm.h - argmax( model.eta_array() )

    #guessed_pssm = model.expected_pssm()
    #assert len( guessed_pssm ) == len( sample.pssm.theta )
    #KL = sum(
    #       [
    #               discrete_KL( t, g )
    #               for t, g
    #               in zip( sample.pssm.theta[1:], guessed_pssm[1:] )
    #       ]
    #)
    #print 'KL per base:', KL / sample.pssm.K

def generate_test(
        N,
        K,
        av_length,
        seed = 1
):
    import time
    start = time.clock()
    print 'Seeding:', seed
    numpy.random.seed( seed )
    phi = 0.05 * numpy.ones( 4, dtype = numpy.float64 )
    varphi = 100.0 * numpy.ones( 4, dtype = numpy.float64 )
    alpha = [ 7.0, 7.0 ] # try and make sure we have even spread of gap and non-gap data
    sample = ModelSample(
            N = N,
            K = K,
            av_length = av_length,
            phi = phi,
            varphi = varphi,
            alpha = alpha,
            verbose = True
    )
    # print sample
    print 'True gap at position:', sample.pssm.h
    #format_weblogo_from_dist(
    #       sample.pssm.theta[1:,:],
    #       'pssms/%s_real' % seed,
    #       'png',
    #)
    dist_of_sites = sample.dist_of_sites()
    format_weblogo_from_dist(
            sample.dist_of_sites(),
            'pssms/%s_data' % seed,
            'png',
    )
    model = create_model(
            K,
            [
                    ''.join( index_to_base( i ) for i in s )
                    for s in sample.seqs
            ]
    )
    learn_model( model, 'pssms/%s_learnt' % seed )
    print 'Took %fs' % (time.clock() - start)
    compare_model_with_truth( model, sample )

def test_with_synthetic():
    for test_id in xrange(1, 20):
        generate_test(
                N = 20,
                K = 10,
                av_length = 40,
                seed = test_id
        )

def test_C_impl():
    model = VariationalModel_C.create(
            4,
            [
                    'aaggt',
                    'aaggc',
                    'aaggg',
                    'tcctt',
                    'gcctt',
                    'ccctt',
                    # 'ttttcg',
                    # 'aaaagc',
                    # 'aaaatt',
                    # 'aagaac',
                    # 'taataa',
            ]
    )
    learn_model( model, 'pssms/C_impl_learnt' )
    return model

def test_with_chip_chip(
        fragment = 'T00594',
        num_pssms = 1,
        max_updates = 2,
        K = 10
):
    print 'Fragment = %s' % fragment
    import corebio.seq_io.fasta_io
    seqs = [
            str( s )
            for s
            in corebio.seq_io.fasta_io.iterseq(
                    open( 'fragments/%s.fa.masked' % fragment, 'r' ),
                    corebio.seq.dna_alphabet
            )
    ]
    # seqs = seqs[:20]
    print '%d sequences' % len(seqs)
    print '%d bases' % sum( [ len(s) for s in seqs ] )
    # model = VariationalModel(
    #       [ seq_to_numeric( seq ) for seq in seqs[:20] ],
    #       4,
    #       alpha = [ 1.0 ] * 2,
    #       phi = [ 0.1 ] * 4,
    #       varphi = [ 10.0 ] * 4,
    # )
    model = VariationalModel_C.create(
            K,
            seqs
    )
    for i in xrange( num_pssms ):
        print 'Learning pssm: %d' % i
        learn_model(
                model,
                'pssms/%s_%d' % (fragment, i),
                max_updates = max_updates
        )
        blanked = model.blank_sites( 0.5 )
        print 'Blanked %d sites' % blanked
        model.initialise_variational_params()
    return model

def test_all_chip_chip():
    # ordered by size
    for fragment in [
            'T00594',
            'T00163',
            'T00759',
            'T00368',
            'T03286',
            'T00140',
            'T00671',
            'T09363',
            'T03828',
            'T08969',
            'T00781'
    ]:
        import time
        start = time.clock()
        test_with_chip_chip( fragment = fragment )
        seconds = time.clock() - start
        print 'Took %ds (=%d minutes)' % ( int(seconds), int(seconds/60) )

def test_compressed_dna_performance():
    setup = 'from biopsy.gapped_pssms import TestCompressedDnaSeq; test = TestCompressedDnaSeq( 10000000 )'
    number = 10
    from timeit import Timer
    t = Timer( "test.test_char()", setup )
    print 'Char:', t.timeit( number = number )
    t = Timer( "test.test_unsigned()", setup )
    print 'Unsigned:', t.timeit( number = number )
    t = Timer( "test.test_compressed()", setup )
    print 'Compressed:', t.timeit( number = number )

if '__main__' == __name__:
    import time
    start = time.clock()
    # test_compressed_dna_performance()
    # test_all_chip_chip()
    # model = test_with_chip_chip()
    model = simple_test()
    # model = test_C_impl()
    # test_with_synthetic()
    seconds = time.clock() - start
    print 'Overall took %ds (=%d minutes)' % ( int(seconds), int(seconds/60) )
