#
# Copyright John Reid 2006
#

from _biopsy import *

def get_constant_score( h ):
    return 1.0



def _dist_str( self ):
    return \
            '%(0)0.7f, %(1)0.7f, %(2)0.7f, %(3)0.7f' % {
                    '0' : self.get( 0 ),
                    '1' : self.get( 1 ),
                    '2' : self.get( 2 ),
                    '3' : self.get( 3 )
            }
NucleoDist.__str__ = _dist_str

def _pssm_str( self ):
    return "\n".join( str( d ) for d in self )
Pssm.__str__ = _pssm_str





def fill_pssm_cache_from_transfac():
    """Fills pssm cache with transfac values"""
    import transfac
    print 'Pseudo counts:', PssmParameters.singleton().pseudo_counts
    for acc in get_transfac_pssm_accessions( transfac.PssmFilter.all_pssms() ):
        get_pssm( acc )
        # save_pssm_cache_state()
    save_pssm_cache_state()


def matrix_length_histogram( pssm_accs = None ):
    """
    Displays a histogram of matrix lengths
    """
    from pylab import hist, show, figure, title
    if not pssm_accs:
        import transfac
        # take all transfac pssms by default
        pssm_accs = get_transfac_pssm_accessions( transfac.PssmFilter.all_pssms() )
    for acc in pssm_accs:
        if len( get_pssm( acc ).pssm ) > 40: print acc
    lengths = [ len( get_pssm( acc ).pssm ) for acc in pssm_accs ]
    figure()
    hist( lengths, bins = max( lengths ), width = 1, normed = True )
    title( 'Transfac matrix lengths' )
    show()

def matrix_num_seqs_histogram( pssm_accs = None ):
    """
    Displays a histogram of then number of sequences matrices are built from
    """
    from pylab import hist, show, figure, title
    if not pssm_accs:
        # take all transfac pssms by default
        pssm_accs = get_transfac_pssm_accessions( PssmFilter.all_pssms() )
    lengths = []
    for acc in pssm_accs:
        l = len( get_transfac_pssm_sequences( acc ) )
        if l > 0:
            lengths.append( l )
        if l < 0:
            print acc
    figure()
    hist( lengths, bins = max( lengths ), width = 1, normed = True )
    title( '# sequences matrices built from' )
    show()






def get_pssm_likelihoods( p ):
    """
    Returns likelihoods, odds ratios and p(binding) lists for pssm.

    Returns a tuple:
    (
            bind,
            back,
            cum_bind,
            cum_back,
            odds_ratio,
            cum_odds_ratio,
            p_bind,
            cum_p_bind,
            p_value_p_bind
    )
    """
    bind = p.get_dist( True, False )
    back = p.get_dist( False, False )
    cum_bind = p.get_dist( True, True )
    cum_back = p.get_dist( False, True )
    odds_ratio = []
    cum_odds_ratio = []
    p_bind = []
    cum_p_bind = []
    p_value_p_bind = []
    for i in range( len( bind ) ):
        score = float( i ) / float( len( bind ) - 1 )
        odds_ratio.append(
                get_odds_ratio(
                        score,
                        bind,
                        back ) )
        cum_odds_ratio.append(
                get_odds_ratio(
                        score,
                        cum_bind,
                        cum_back ) )
        p_bind.append( get_p_binding( odds_ratio[i] ) )
        cum_p_bind.append( get_p_binding( cum_odds_ratio[i] ) )
        p_value_p_bind.append( get_p_binding_using_p_value( score, cum_bind ) )
    return (
            bind,
            back,
            cum_bind,
            cum_back,
            odds_ratio,
            cum_odds_ratio,
            p_bind,
            cum_p_bind,
            p_value_p_bind
    )

def get_pssm_likelihoods_for_score( p, score ):
    """
    Returns likelihoods, odds ratios and p(binding) lists for particular pssm
    score.

    Returns a tuple:
    (
            bind,
            back,
            cum_bind,
            cum_back,
            odds_ratio,
            cum_odds_ratio,
            p_bind,
            cum_p_bind,
            p_value_p_bind
    )
    """
    (
            bind,
            back,
            cum_bind,
            cum_back,
            odds_ratio,
            cum_odds_ratio,
            p_bind,
            cum_p_bind,
            p_value_p_bind
    ) = get_pssm_likelihoods( p )
    idx = get_likelihood_index( len( bind ), score )
    return (
            bind[ idx ],
            back[ idx ],
            cum_bind[ idx ],
            cum_back[ idx ],
            odds_ratio[ idx ],
            cum_odds_ratio[ idx ],
            p_bind[ idx ],
            cum_p_bind[ idx ],
            p_value_p_bind[ idx ]
    )



def create_cumulative_dists( dist ):
    cum = 1.0
    result = Likelihoods()
    for i in range( len( dist ) ):
        result.append( cum )
        cum -= dist[ i ]
    return result


def plot_likelihoods( p, name, score = None ):
    from pylab import plot, show, cla, clf, legend, figure, xlabel, ylabel, \
            title, axvspan
    old_odds = PssmParameters.singleton().binding_background_odds_prior
    PssmParameters.singleton().binding_background_odds_prior = 1.0
    (
            bind,
            back,
            cum_bind,
            cum_back,
            odds_ratio,
            cum_odds_ratio,
            p_bind,
            cum_p_bind,
            p_value_p_bind
    ) = get_pssm_likelihoods( p )
    scores = [ float(i) / (len(bind) - 1.0) for i in range( len( bind ) ) ]
    # cla()
    # clf()
    figure()
    plot( scores, bind, 'g-', label='binding' )
    plot( scores, back, 'b-', label='background' )
    plot( scores, cum_bind, 'g--', label='binding (cumulative)' )
    plot( scores, cum_back, 'b--', label='background (cumulative)' )
    plot( scores, p_bind, 'r-', label='p(binding)' )
    plot( scores, cum_p_bind, 'r--', label='p(binding) (cumulative)' )
    plot( scores, p_value_p_bind, 'y--', label='p(binding) (p-value)' )
    # legend( loc='center left' )
    xlabel( 'score' )
    title( name )
    if score:
        idx = get_likelihood_index( len( bind ), score )
        axvspan( score, score )
    show()
    PssmParameters.singleton().binding_background_odds_prior = old_odds
