#
# Copyright John Reid 2006
#

import cPickle, os, math, scipy.stats.distributions
from binding_hit import *
from families import *
from graph import *
from paralogs import *
from remo import *

class Bunch( object ):
    """Holds a bunch of objects as attributes"""

    def __init__(self, **kwds):
        """Takes named arguments, one for each object"""
        self.__dict__.update(kwds)

    def __str__(self):
        return '\n'.join(
                [
                        "%s=%r" % (attribute, value)
                        for (attribute, value)
                        in self.__dict__.iteritems()
                ]
        )


def make_binder_pair( h1, h2 ):
    """Returns a binder pair from the two hits

    A pair has the form:

            ( binder1, binder2, orientation1, orientation2 )

    Where the binders are pssm names typically and the orientation describes
    the strand the pssm binds to (true for positive strand).
    """
    b1 = h1.binder
    b2 = h2.binder
    s1 = h1.location.positive_strand
    s2 = h2.location.positive_strand

    # first make sure binders are in correct order
    if b1 > b2:
        # they aren't so swap them
        b1, b2 = b2, b1
        s1, s2 = not s2, not s1

    return (
            b1,
            b2,
            s1,
            s2
    )

def pair_analysis(
        remome,
        analysis,
        max_distance,
        hit_filter = None,
        pair_filter = None,
        distance_filter = None,
        location_filter = None
):
    """Examines the analysis for pairs satisfying the distance and filter reqs.

    max_distance is the maximum separation that is looked for.
    The filters are called to check if the hits, pairs and or location are of
    interest. They are disregarded if not.
    """
    if isinstance( remome, str ):
        print 'Loading remome from:', remome
        remome = Remome.deserialise( remome )

    if isinstance( analysis, str ):
        print 'Loading analysis from:', analysis
        analysis = Analysis.deserialise( analysis )

    # maps binder pairs to vectors of counts indexed by max_distance
    pair_dists = { }

    num_analyses = 0
    num_hits = 0
    num_pairs = 0

    # for each piece of analysis
    for k in analysis.get_keys():
        # print k
        try:
            remo, alignment = Remome.get_remo_from_id( remome, k )
        except RuntimeError, e:
            # print "Could not get remo:", e
            continue

        seq_len = len( remo.get_sequence_for( alignment.centre_sequence, True ) )
        num_analyses += 1
        hits = analysis.get_hits_for( k )
        for h1 in hits:
            if None != hit_filter and not hit_filter( h1 ): continue
            num_hits += 1

            for h2 in hits:
                if None != hit_filter and not hit_filter( h2 ): continue
                # how far apart are they?
                distance = h2.location.start() - h1.location.end()

                # is h1 before h2 with no overlap?
                if distance < 0: continue # no

                # are they within the max_distance?
                if distance >= max_distance: continue # no

                # is the first one too close to the end?
                if seq_len - ( h1.location.end() + h2.location.length ) < max_distance: continue # yes

                # is the second one too close to the beginning?
                if h2.location.start() - h1.location.length < max_distance: continue # yes

                if None != distance_filter and not distance_filter( distance ): continue

                binder_pair = make_binder_pair( h1, h2 )

                # check that we want to consider this pair in this remo
                if None != pair_filter and not pair_filter( binder_pair ): continue
                num_pairs += 1

                # check that we want to consider this location in this remo
                if (
                        None != location_filter
                        and not location_filter( binder_pair, remo, alignment )
                ):
                    continue

                # make sure we have an entry in the map
                if not pair_dists.has_key( binder_pair ):
                    pair_dists[ binder_pair ] = 0
                    pair_dists[ binder_pair ] = [ 0 ] * max_distance
                dist = pair_dists[ binder_pair ]

                # increase the count for this distance
                dist[ distance ] += 1


    print 'Looked at %d analyses' % ( num_analyses )
    print 'Looked at %d hits' % ( num_hits )
    print 'Looked at %d different combinations of factors' % ( len( pair_dists ) )
    print 'Looked at %d instances of these pairs' % ( num_pairs )

    return pair_dists

class GraphFilter( object ):
    """Checks if putative interesting hits have already been seen in other
    genes that are connected in the supplied graph.

    Designed to be used as a location filter in pair_analysis()
    """
    def __init__( self, g ):
        # a graph of family relationships
        self.g = g
        # we need a mapping from pairs to genes in which the pair has already been
        # seen
        self.genes_for_pair = { }

    def __call__( self, binder_pair, remo, alignment ):
        # get the gene id
        gene = str( alignment.centre_sequence ).split()[0]

        # does our graph know about this gene?
        if not graph_has_vertex( self.g, gene ):
            return False
        v = graph_vertex( self.g, gene )

        # get the set of genes for which we have already seen this pair
        if binder_pair not in self.genes_for_pair:
            gene_set = set()
            self.genes_for_pair[ binder_pair ] = gene_set
        else:
            gene_set = self.genes_for_pair[ binder_pair ]

        # are any of them paralogs of this gene?
        for g2 in gene_set:
            if graph_are_connected( self.g, v, graph_vertex( self.g, g2 ) ):
                return False
        # otherwise - add to gene set and return true
        gene_set.add( gene )
        return True

class AndFilter( object ):
    """Combines two other filters with a logical AND"""
    def __init__( self, filter1, filter2 ):
        self.filter1 = filter1
        self.filter2 = filter2

    def __call__( self, binder_pair, remo, alignment ):
        return self.filter1(
                binder_pair, remo, alignment
        ) and self.filter2(
                binder_pair, remo, alignment
        )

class OrFilter( object ):
    """Combines two other filters with a logical OR"""
    def __init__( self, filter1, filter2 ):
        self.filter1 = filter1
        self.filter2 = filter2

    def __call__( self, binder_pair, remo, alignment ):
        return self.filter1(
                binder_pair, remo, alignment
        ) or self.filter2(
                binder_pair, remo, alignment
        )

def get_default_pair_filter(
        family_file = None,
        paralog_file = None,
        subtypes = None
):
    """Get the default filter for filtering pairs from hits

    family_file is a file that defines protein families
    paralog_file is a file that defines paralogs
    subtypes is an array of paralog types that are used

    Will ignore hits that are match other hits from the same families or
    from genes that are paralogs.
    """
    if None == family_file: family_file = 'C:/Data/Ensembl/mouse_families.txt'
    if None == paralog_file: paralog_file = 'C:/Data/Ensembl/mouse_paralogs.txt'
    if None == subtypes:
        subtypes = [
                # 'Fungi/Metazoa group',
                # 'Coelomata',
                # 'Bilateria',
                # 'Chordata',
                # 'Amniota',
                # 'Tetrapoda',
                'Euteleostomi',
                'Theria',
                'Eutheria',
                'Euarchontoglires',
                'Murinae',
                'Mus musculus',
        ]

    return AndFilter(
            GraphFilter(
                    graph_generate(
                            read_families_file(
                                    family_file
                            )
                    )
            ),
            GraphFilter(
                    graph_generate(
                            read_paralog_file(
                                    paralog_file,
                                    subtypes
                            )
                    )
            )
    )


def calc_and_save_pair_analysis(
        remome_file,
        analysis_file,
        max_distance,
        filename
):
    pair_dists = pair_analysis(remome_file, analysis_file, max_distance)
    save_pair_analysis( pair_dists, filename )



def log_factorial( n ):
    if 0 > n: raise RuntimeError, 'n must be non-negative: %f' % ( n )
    if 0 == n: return 0
    return math.log( n ) + log_factorial( n - 1 )

def log_n_choose( y ):
    n = sum( y )
    return (
            log_factorial( n )
            - sum( log_factorial( x ) for x in y )
    )

def log_gamma( x ):
    import scipy.special
    try:
        return math.log( scipy.special.gamma( x ) )
    except:
        print 'x: %f' % ( x )
        print 'scipy.special.gamma( x ): %f' % ( scipy.special.gamma( x ) )
        raise

def log_Z( alpha, m ):
    return (
            sum( log_gamma( alpha * m_i ) for m_i in m )
            - log_gamma( alpha )
    )

_gen = scipy.stats.distributions.gamma_gen( name = 'gamma' )
def dirichlet_rv( alpha, m ):
    y = [ _gen.rvs( alpha * m_i )[0] for m_i in m ]
    s = sum( y )
    return [ y_i / s for y_i in y ]




class MultinomialSample:
    """A sample from a multinomial distribution. Compares models.
    """

    def __init__( self, y ):
        for y1 in y: assert isinstance( y1, int )
        self.y = y
        self.n = sum( y )
        self.k = len( y )
        self.log_n_choose = log_n_choose( y )

    def ll_under_multinomial( self, p = None ):
        if None == p: p = [ 1.0 / self.k ] * self.k
        return (
                self.log_n_choose
                + sum( y_i * math.log( p_i ) for p_i, y_i in zip( p, self.y ) )
        )

    if False: # old implementation with MathRangeErrors
        def ll_under_dirichlet( self, alpha = 1.0, m = None ):
            if None == m: m = [ 1.0 / self.k ] * self.k
            a = [ alpha * m_i + y_i for m_i, y_i in zip( m, self.y ) ]
            alpha_prime = sum( a )
            m_prime = [ a_i / alpha_prime for a_i in a ]
            return (
                    self.log_n_choose
                    - log_Z( alpha, m )
                    + log_Z( alpha_prime, m_prime )
            )
    else:
        def ll_under_dirichlet( self, alpha = 1.0, m = None ):
            if None == m: m = [ 1.0 / self.k ] * self.k
            a = [ alpha * m_i for m_i in m ]
            return (
                    self.log_n_choose
                    - sum(
                            math.log( alpha + j )
                            for j
                            in range( self.n )
                    )
                    + sum(
                            sum(
                                    math.log( a[i] + j )
                                    for j in range( self.y[i] )
                            )
                            for i in range( self.k )
                    )
            )

    def model_comparison( self, p = None, alpha = 1.0, m = None ):
        """Equivalent of ll_under_multinomial() - ll_under_dirichlet()"""
        if None == p: p = [ 1.0 / self.k ] * self.k
        if None == m: m = [ 1.0 / self.k ] * self.k
        a = [ alpha * m_i for m_i in m ]
        sum_a = sum( a )
        return (
                sum(
                        y_i * math.log( p_i )
                        for p_i, y_i
                        in zip( p, self.y )
                )
                + sum(
                        math.log( alpha + j )
                        for j
                        in xrange( self.n )
                )
                - sum(
                        sum(
                                math.log( a[i] + j )
                                for j
                                in xrange( self.y[i] )
                        )
                        for i
                        in xrange( self.k )
                )
        )

def display_pair_distribution( p ):
    """Shows a histogram of the separation distribution for a pair"""
    from pylab import figure, bar, show, title, xlabel, ylabel
    figure()
    bar( range( len( p.dist ) ), p.dist )
    title(
            '%s : %s - %s\n# samples = %d, LOR = %f' % (
                    str( p.binder_pair ),
                    get_transfac_pssm_name( p.binder_pair[0] ),
                    get_transfac_pssm_name( p.binder_pair[1] ),
                    sum( p.dist ),
                    p.log_odds_ratio,
            )
    )
    xlabel( 'separation' )
    ylabel( 'count' )
    # show()

def sort_pairs_by_log_odds( pairs ):
    pairs.sort( key = lambda p: p.log_odds_ratio )
    return pairs

def calculate_pair_log_odds( pairs_dists ):
    """Calculates the log odds for a pair

    Returns a sequence of Bunch objects containing the binder_pair, the
    log_odds_ratio and the distribution of distances
    """
    analyses = [ ]
    for binder_pair, dist in pairs_dists.iteritems():
        analyses.append(
                Bunch(
                        binder_pair = binder_pair,
                        log_odds_ratio = MultinomialSample(
                                dist
                        ).model_comparison( alpha = float( len( dist ) ) ),
                        dist = dist
                )
        )

    # sort by log odds
    analyses.sort( key = lambda p: p.log_odds_ratio )
    return analyses

def write_pair_separation_histograms( pairs, directory ):
    """Writes separarion histograms for the pairs into the given directory

    May not work if matplotlib or pylab are already imported as it needs to set
    use options for matplotlib before pylab is imported
    """
    import matplotlib
    matplotlib.use('Agg')
    import pylab, biopsy
    if not os.access( directory, os.W_OK ): os.mkdir( directory )
    for i, p in enumerate( pairs ):
        display_pair_distribution( p )
        pylab.savefig(
                os.path.join(
                        directory,
                        '%03d' % ( i )
                )
        )

def get_pairs_for_pssm( pairs, pssm ):
    """Returns a sequence of pairs which the given pssm is in
    """
    return [
            p
            for p
            in pairs
            if (
                    p.binder_pair[0] == pssm
            ) or (
                    p.binder_pair[1] == pssm
            )
    ]


def get_gene_universe( genes_file ):
    """Builds a unique sorted list of genes from gene file"""
    if isinstance( genes_file, str ):
        genes_file = open( genes_file, 'r' )
    return set( [ l.strip() for l in genes_file ] )


def print_pairs_go_analysis( pairs, gene_universe, analysis ):
    """Prints the go analysis for the given pairs
    """
    import rpy, r_go

    mart = r_go.get_mart( "mmusculus_gene_ensembl" )
    categs = r_go.categorise_genes( mart, gene_universe )

    #pair = ( 'M00349', 'M00350', True, False )
    for p in pairs:
        pair = p.binder_pair
        print pair

        in_analysis = find_pair_in_analysis(
                analysis,
                pair,
                max_separation = 45
        )

        genes = [
                seq.split(' ')[0]
                for seq, hits
                in in_analysis.iteritems()
        ]

        print len(genes)
        print "\n".join( genes )

        result = rpy.with_mode(
                rpy.NO_CONVERSION,
                rpy.r.analyseGoAnnotations
        )(
                categs,
                genes
        )

        rpy.r.printAnnotationResult( result, 10 )

def get_pairs_from( all_pairs, pssms_1, pssms_2 ):
    """Returns pairs from all_pairs that have one pssm from pssms_1 and one pssm
    from pssms_2
    """
    return [
            p
            for p
            in all_pairs
            if (
                    (p.binder_pair[0] in pssms_1) and (p.binder_pair[1] in pssms_2)
            ) or (
                    (p.binder_pair[1] in pssms_1) and (p.binder_pair[0] in pssms_2)
            )
    ]

def research_pair(
        analysis,
        categs,
        pair,
        max_separation = 45
):
    """Look for infomation about a particular pair

    analysis: The analysis file or object
    categs: Go categories R object
    pair: The definition of the pair
    """
    import rpy

    in_analysis = find_pair_in_analysis(
            analysis,
            pair,
            max_separation
    )
    genes = [ s.split()[0] for s in in_analysis ]
    print '%s is in:\n%s' % ( str(pair), "\n".join( genes ) )

    result = rpy.with_mode(
            rpy.NO_CONVERSION,
            rpy.r.analyseGoAnnotations
    )(
            categs,
            genes
    )
