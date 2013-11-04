
try:
    import _hmm as hmm
except:
    print 'Could not import _hmm, trying hmm'
    import hmm
import numpy, os, math


nt = hmm.numpy_test()
nt.print_shape()
nt.array = numpy.zeros( (4,10) )
nt.print_shape()
nt.print_first()
nt.array = 4.5 * numpy.ones( (4,10) )
nt.print_first()

print nt.array

raise RuntimeError()

def create_pssm_model():
    position_distributions = numpy.array(
            [
                    hmm.dirichlet_draw( 0.01 * numpy.ones( 4 ) )
                    for i in xrange( K )
            ]
    )
    return hmm.build_pssm_model(
            p_binding_site = p_binding_site,
            position_distributions = numpy.array(
                    [
                            hmm.dirichlet_draw( 10.0 * numpy.ones( 4 ) )
                            for i in xrange( K )
                    ]
            )
    )

def create_fully_connected_model():
    return hmm.build_fully_connected_hmm( K, 4 )

def generate_sequences( generating_model ):
    seqs = [ ]
    true_states = [ ]
    for i in xrange( N ):
        states, seq = generating_model.sample( T )
        seq.dtype = numpy.int32
        seqs.append( seq )
        true_states.append( states )
    return seqs, true_states

def do_baum_welch( model, seqs ):
    def bw_callback( LL ):
        #print LL
        for i in xrange( N ):
            #print '%d: %f'  % ( i, model.forward( seqs[ i ] )[ 0 ] )
            pass
        #print model.B
    #print model.B
    model.baum_welch(
            seqs,
            max_iterations = 10,
            callback = bw_callback
    )

K = 2
p_binding_site = 0.1
for N in [ 1, 2, 10 ]:
    for T in [ 10, 1000 ]:
        for model_generator in create_fully_connected_model, create_pssm_model:
            print N, T, model_generator
            for seed in xrange( 1, 101 ):
                #print seed
                seqs, true_states = generate_sequences( create_pssm_model() )
                model = model_generator()
                do_baum_welch( model, seqs )
