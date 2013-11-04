#
# Copyright John Reid 2006
#

import biopsy, sys

# filename = sys.argv[1]

def load_remome():
    filename = 'c:/data/remos/100/100.cleaned'
    print 'Loading remome from %s' % ( filename )
    return biopsy.Remome.deserialise( filename )


def pair_log_odds_histogram():
    pairs = biopsy.load(
            'c:/data/remos/100/100.pairs_lor',
            display = 'pairs log odds'
    )
    import pylab, math
    pylab.hist(
            [ p.log_odds_ratio for p in pairs ],
            100
    )
    pylab.show()
# pair_log_odds_histogram()

def test_get_pairs_for_pssm():
    pssm = 'M00736'
    pairs = biopsy.load(
            'c:/data/remos/100/100.top_pairs',
            display = 'pairs log odds'
    )
    pssm_pairs = biopsy.get_pairs_for_pssm(
            pairs,
            pssm
    )
    print 'Got %d pairs for %s from total of %d' % \
            ( len(pssm_pairs), pssm, len(pairs) )
    for p in pssm_pairs[:20]:
        print p.binder_pair, p.log_odds_ratio
    biopsy.write_pair_separation_histograms( pssm_pairs[:20], pssm )
# test_get_pairs_for_pssm()

def test_pair_go_analysis():
    import rpy

    rpy.r.library("Category")
    rpy.r.library("biomaRt")
    rpy.r.source('../../R/go_categorise.R')

    remome_threshold = 100
    analysis_file = 'c:/data/remos/%d/%d.analysis' \
            % ( remome_threshold, remome_threshold )
    pairs_file = 'C:/Data/ReMos/%d/%d.top_pairs' \
            % ( remome_threshold, remome_threshold )
    gene_universe_file = 'c:/Data/remos/%d/%d.uniq_genes' \
            % ( remome_threshold, remome_threshold )

    gene_universe = rpy.r.read_table( gene_universe_file )["V1"]
    analysis = biopsy.Analysis.deserialise( analysis_file )
    pairs = biopsy.load( pairs_file, display = 'pairs' )

    biopsy.print_pairs_go_analysis(
            pairs[:150],
            gene_universe,
            analysis
    )
    print 'done'
# test_pair_go_analysis()

def look_for_factor_pair( name, all_pairs, pssms1, pssms2 ):
    both_pairs = biopsy.get_pairs_from( all_pairs, pssms1, pssms2 )
    biopsy.sort_pairs_by_log_odds( both_pairs )
    print 'Found %d pairs' % len( both_pairs )
    for p in both_pairs:
        print p.binder_pair, p.log_odds_ratio
    biopsy.write_pair_separation_histograms( both_pairs, name )
def test_look_for_factor_pair():
    print 'test_look_for_factor_pair'
    all_pairs = biopsy.load(
            'c:/data/remos/0/0.top_pairs',
            display = 'pairs log odds'
    )
    sascha_factors = [
            "T00140", # c-Myc                       human, Homo sapiens
            "T00163", # CREB                        human, Homo sapiens
            "T00368", # HNF-1alpha-A        human, Homo sapiens
            "T00594", # RelA-p65            human, Homo sapiens
            "T00671", # p53                         human, Homo sapiens
            "T00759", # Sp1                         human, Homo sapiens
            "T00781", # TAF(II)250          human, Homo sapiens
            "T03286", # HNF-6alpha          human, Homo sapiens
            "T03828", # HNF-4alpha          human, Homo sapiens
    ]
    sascha_pssms = set()
    for factor in sascha_factors:
        pssms_for_factor = biopsy.get_pssms_for_factor( factor )
        print 'Pssms for %s: %s' % ( factor, " ".join( pssms_for_factor ) )
        for pssm in pssms_for_factor:
            sascha_pssms.add( pssm )
    look_for_factor_pair(
            'sascha',
            all_pairs,
            sascha_pssms,
            sascha_pssms
    )
test_look_for_factor_pair()

def test_find_pair_in_analysis():
    remome_threshold = 0
    analysis_file = 'c:/data/remos/%d/%d.analysis' % ( remome_threshold, remome_threshold )
    analysis = biopsy.Analysis.deserialise( analysis_file )
    hits = biopsy.find_pair_in_analysis(
            analysis,
            ( 'M01019', 'M01043', True, True ),
            max_separation = 45
    )
    print 'Got %d hits' % len( hits )
    for s, v in hits.iteritems():
        print s
        for hit in v:
            print max(
                    hit[1].location.start() - hit[0].location.end(),
                    hit[0].location.start() - hit[1].location.end()
            )
# test_find_pair_in_analysis()

def test_pair_analysis():
    import webbrowser
    remome_threshold = 100
    remome_file = 'c:/data/remos/%d/%d.cleaned' % ( remome_threshold, remome_threshold )
    analysis_file = 'c:/data/remos/%d/%d.analysis' % ( remome_threshold, remome_threshold )
    location_filter = biopsy.get_default_pair_filter()
    hits = biopsy.pair_analysis(
            remome_file,
            analysis_file,
            max_distance = 45,
            distance_filter = lambda d: d == 25,
            pair_filter = lambda p: p == ( 'M00471', 'M00791', False, False ),
            location_filter = location_filter
    )
    for k, v in hits.iteritems():
        print k
        print v
        for g in location_filter.filter2.genes_for_pair[ k ]:
            print g
            webbrowser.open( biopsy.ensembl_url( g ) )
        print
# test_pair_analysis()

def test_research_pair():
    import biopsy.r_go, rpy

    pair = ('M00293', 'R04602', False, True)
    remome_threshold = 100
    analysis_file = 'c:/data/remos/%d/%d.analysis' % ( remome_threshold, remome_threshold )
    pairs_file = 'C:/Data/ReMos/%d/%d.top_pairs' % ( remome_threshold, remome_threshold )
    gene_universe_file = 'c:/Data/remos/%d/%d.uniq_genes' % ( remome_threshold, remome_threshold )

    gene_universe = rpy.r.read_table( gene_universe_file )["V1"]
    analysis = biopsy.Analysis.deserialise( analysis_file )
    pairs = biopsy.load( pairs_file, display = 'pairs' )

    mart = biopsy.r_go.get_mart( "mmusculus_gene_ensembl" )
    categs = biopsy.r_go.categorise_genes( mart, gene_universe )
    biopsy.research_pair( analysis, mart, categs, pair )
# test_research_pair()

def saschas_pair_filter( binder_pair ):
    return (
            ( binder_pair[0]=='M01019' or binder_pair[1]=='M01043' )
            or ( binder_pair[1]=='M01019' or binder_pair[0]=='M01043' )
    )

def find_saschas_pairs():
    pair_dists = biopsy.pair_analysis(
            'c:/data/remos/0/0.cleaned',
            'c:/data/remos/0/0.analysis',
            45,
            pair_filter = saschas_pair_filter,
            location_filter = biopsy.get_default_pair_filter(
                    family_file = 'C:/Data/Ensembl/mouse_families.txt',
                    paralog_file = 'C:/Data/Ensembl/mouse_paralogs.txt',
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
            )
    )
    biopsy.save( pair_dists, 'saschas.pairs' )
# find_saschas_pairs()

def test_saschas_pairs():
    pair_dists = biopsy.load( 'saschas.pairs' )
    for p, d in pair_dists.iteritems():
        if ( p[0] == 'M01019' and p[1] == 'M01043' ) \
                or ( p[1] == 'M01019' and p[0] == 'M01043' ):
            print p
            print d
# test_saschas_pairs()

# test find sequence in remome
def test_find_sequence_in_remome():
    for aligned, remo in biopsy.Remome.deserialise(
            'c:/data/remos/0/0.cleaned'
    ).find_sequence(
            # 'AGTCCAAGGTGTTGAATGTTGCCACTTCAAGCCTAAACTTTCTAGGAACACCTAAGTGGGTGGCAGCTT',
            'AGTCCAAG',
            masked = False,
            examine_all_species = True
    ):
        print aligned.get_remo_id( remo )
# test_find_sequence_in_remome()
