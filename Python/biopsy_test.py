#
# Copyright John Reid 2006
#

import sys, time, profile, math, random, pickle, unittest, biopsy
print biopsy.get_build()
# print dir( biopsy )


def test_():
    print '**************** test_ ****************'
# test_()



def test_remome_analysis():
    print '**************** test_remome_analysis ****************'
    analysis = biopsy.RemomeAnalysis( biopsy.Remome.load('C:/Data/ReMos/remo_space.bin') )
    analysis.analysis.serialise( 'analysis.bin' )
    analysis_copy = biopsy.Analysis.deserialise( 'analysis.bin' )
    pssms = biopsy.get_transfac_pssm_accessions( biopsy.get_default_transfac_pssm_filter() )
    threshold = 0.05
    phylo_threshold = 0.02
    try:
        analysis.analyse(
                pssms,
                threshold,
                phylo_threshold
        )
        analysis.analysis.serialise( 'analysis.bin' )
    except:
        analysis.analysis.serialise( 'analysis.bin' )
        raise

test_remome_analysis()



def test_remome():
    print '**************** test_remome ****************'
    remome = biopsy.Remome.load('C:/Data/ReMos/remo_space.bin')
    seqs = biopsy.Remome.make_gene_alignment_map(remome.get_aligned_sequences())
    aligned = remome.get_aligned_sequences()
    gene_alignment = biopsy.Remome.make_gene_alignment_map(aligned)
    print len( gene_alignment.get_genes( 'ENSMUSG' ) )
    seqs = gene_alignment.get_alignments_for(
            biopsy.Remome.EnsemblId.parse( 'ENSMUSG00000035799' ) )
    print len(seqs)
    for s in seqs:
        print s
        for r in remome.get_remos_for(s):
            print biopsy.Remome.get_remo_id( s, r )
            # print '\t', s.centre_sequence, '\n\t', '\n\t'.join( str(i) for i in r.get_sequence_ids() )
            print '\t', len( r.get_sequence_ids() ),  r.get_sequence_for( s.centre_sequence, True )
            break
        print
        # break
# test_remome()

def test_pssm_pseudo_counts():
    sascha_pssms = biopsy.SequenceVec()
    sascha_acc = 'M00975'
    # sascha_seq = 'gtaaaccaggctgcctGAgaacttgttgcgaatcc'
    sascha_seq = 'ttgttgcga'
    sascha_seq = 'ttgttgcaa'
    # plot_likelihoods( biopsy.get_pssm( 'M00975' ), 'M00975' )
    # plot_likelihoods( biopsy.get_pssm( 'R02146' ), 'R02146' )
    print 'Binding,Background,odds,p(binding),cumulative p(binding),Sequence'
    biopsy.PssmParameters.singleton().use_p_value = True;
    # biopsy.PssmParameters.singleton().binding_background_odds_prior = 1;
    for pc in [ 0.0, 0.25, 0.5, 1.0, 2.0 ]:
        # force cache load
        biopsy.get_pssm( sascha_acc )
        biopsy.clear_pssm_cache()
        biopsy.PssmParameters.singleton().pseudo_counts = pc
        p = biopsy.get_pssm( sascha_acc )
        score = biopsy.score_pssm( p.pssm, sascha_seq )
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
        ) = biopsy.get_pssm_likelihoods_for_score( p, score )
        print pc,
        print \
                '%f,%f,%f,%f,%f,%f,%f' \
                % \
                ( bind, back, cum_bind, cum_back, p_bind, cum_p_bind, p_value_p_bind )
        biopsy.plot_likelihoods( p, sascha_acc + ': ' + str( pc ), score )
        # print 'Trying with standard distributions'
        # biopsy.PssmParameters.singleton().use_cumulative_dists = False;
        # hits = biopsy.HitVec()
        # biopsy.score_pssm_on_sequence( sascha_acc, sascha_seq, 0.001, hits )
        # print hits
        print 'Trying with cumulative distributions'
        biopsy.PssmParameters.singleton().use_cumulative_dists = True;
        hits = biopsy.HitVec()
        biopsy.score_pssm_on_sequence( sascha_acc, sascha_seq, 0.001, hits )
        print hits
        print
# test_pssm_pseudo_counts()

def test_load_remome():
    # print dir( biopsy )
    biopsy.load_remome( 'C:/Data/ReMos/remo_space.bin' )
# test_load_remome()


def test_likelihoods_indices():
    p = biopsy.get_pssm( 'M00975' )
    dist = p.get_dist( True, False )
    for s in range( len(dist) ):
        score = float(s)/float(len(dist) - 1)
        idx = biopsy.get_likelihood_index( len(dist), score )
        print idx, score
    for score in [ 0.98, 0.99, 1.0 ]:
        print score, biopsy.get_likelihood_index( len(dist), score )
# test_likelihoods_indices()


def examine_pssm_distributions():
    # read the distributions from disk
    f = open( 'pssm_p_binding_dists.txt', 'r' )
    print 'Reading pssm p(binding) distributions from:', f
    dists = pickle.load( f )
    f.close()

    # print dists.keys()

    def display_with_pylab():
        from pylab import plot, subplot, show, cla, clf, LogLocator, MultipleLocator

        def plot_dist( name ):
            # ax = subplot( 111 )
            values = dists[name].values()
            cumulative = [ 0 ]
            for v in values:
                cumulative.append( cumulative[-1] + v )
            plot( cumulative )
            # ax.xaxis.set_major_locator( LogLocator() )

        cla()
        clf()
        plot_dist( 'all' )
        show()
        return
        for i in range( 100 ):
            k = dists.keys()[i+1]
            if 'all' == k: continue
            print k
            plot_dist( k )
        show()
    display_with_pylab()

    def display_with_hippo():
        import hippo
        if not 'app' in dir():
            app = hippo.HDApp()
            canvas = app.canvas()
        def create_data_rep( dist ):
            print hippo.DataRep.names()
            ds = hippo.ListTuple()
            ds.addColumn( 'p(binding)', dist.keys() )
            ds.addColumn( 'count', dist.values() )
            return hippo.DataRep(
                    'XY Plot',
                    ds,
                    [
                            'p(binding)',
                            'count',
                            'nil'
                    ] )

        # show a graph
        # print score_dists[ pssm_acc ]
        xy = hippo.Display ( "XY Plot" )
        xy.addDataRep( create_data_rep( dists[ 'all' ] ) )
        canvas.addDisplay ( xy )
# examine_pssm_distributions()

def test_pssm_distributions():
    pssm_acc = 'M00750'
    transfac_pssms = biopsy.get_transfac_pssm_accessions( biopsy.get_default_transfac_pssm_filter() )
    print 'Got', len( transfac_pssms ), 'pssms'

    # initialise the distributions
    score_dists = { }
    for tp in transfac_pssms:
        score_dists[ tp ] = { }
        for k in range(1,100):
            score_dists[ tp ][ k ] = 0
    score_dists[ 'all' ] = { }
    for k in range(1,100):
        score_dists[ 'all' ][ k ] = 0

    # parse the mouse chromosome, score the pssms and fill in the distributions
    bases = 0
    seq = ''
    start = time.clock()
    for line in open( 'C:/Data/ensembl/chromosomes/Mus_musculus.NCBIM34.dec.dna.chromosome.1.fa', 'r' ):
        if line.startswith( '>' ): continue
        seq += line.strip( '\r\n' ).replace( 'N', '' )
        # Take 1kb at a time
        if len( seq ) >= 1000:
            hits = biopsy.score_pssms_on_sequence( transfac_pssms, seq, 0.0 )
            for h in hits:
                bin = int( 100.0 * h.p_binding )
                if 0 != bin:
                    # print h.binder, bin
                    score_dists[ h.binder ][ bin ] += 1
                    score_dists[ 'all' ][ bin ] += 1
            bases += len( seq )
            seq = ''
            print 'Bases:', bases
        if bases >= 2800000: break
    elapsed = time.clock() - start
    print 'Scored', len( transfac_pssms ), 'pssms on', bases, 'bases in', elapsed, 'seconds'
    print 'Estimate for mouse chromosome 1 (secs):', elapsed * 190000000 / bases
    print 'Estimate for mouse chromosome 1 (days):', elapsed * 190000000 / bases / 60 / 60 / 24
    print 'Estimate for # bases/hour:', bases * 3600 / elapsed

    # remember scores for later
    f = open( 'pssm_p_binding_dists.txt', 'w' )
    print 'Writing pssm p(binding) distributions to:', f
    pickle.dump(score_dists, f)
    f.close()

# test_pssm_distributions()





def test_score_fasta():
    print '******** test_score_fasta()'
    sequences = biopsy.SequenceVec()
    for name, seq in biopsy.parse_fasta( 'c:/analysis/keiths/msx1/enhancerD.fa' ).iteritems():
        print name, ':', len( seq ), 'bases'
        sequences.append( seq )
        if len( sequences) >= 3: break
    phylo_result = biopsy.score_pssms_on_phylo_sequences(
            biopsy.get_transfac_pssm_accessions( biopsy.get_default_transfac_pssm_filter() ),
            sequences,
            0.05 )
    print 'Max Chain:'
    print phylo_result[ 1 ]
    # print biopsy.sort_hits_by_position( phylo_result[ 0 ] )
    print 'Got', len( phylo_result[ 0 ] ), 'hits from', len( sequences[ 0 ] ), 'bases'
# test_score_fasta()




def test_hippo():
    print 'Creating app'
    app = hippo.HDApp()
    print 'Creating canvas'
    canvas = app.canvas()
    print 'Creating data'
    x = []
    for i in range ( 10000 ) :
        x.append ( random.gauss ( 45, 10 ) )
    print 'Displaying'
    hist = Display ( 'Histogram', ( x, ), ('Gaussian', ) )
    print 'Adding display'
    canvas.addDisplay ( hist )
# test_hippo()


def test_fill_pssm_cache():
    biopsy.fill_pssm_cache_from_transfac()
# test_fill_pssm_cache()




def test_transfac_pssms():
    transfac_pssms = biopsy.get_transfac_pssm_accessions( biopsy.get_default_transfac_pssm_filter() )
    for p in transfac_pssms:
        print p, biopsy.get_transfac_pssm_name( p )
    print 'Have', len( transfac_pssms ), 'transfac pssms'
    for acc in [ 'R19099', 'M00418' ]:
        print acc, biopsy.get_transfac_pssm_name( acc )
        biopsy.get_pssm( acc )
        print 'Under pssm'
        for under_pssm in biopsy.get_pssm( acc ).get_dist( True, False ):
            print under_pssm
        print 'Under background'
        for under_background in biopsy.get_pssm( acc ).get_dist( False, False ):
            print under_background
# test_transfac_pssms()



def test_reverse_complement():
    for seq in [
                    '',
                    'a',
                    'c',
                    'g',
                    't',
                    'A',
                    'C',
                    'G',
                    'T',
                    'taca',
                    'tacatcatctgt',
                    'tacatcatctgtctgcagtagtctaacc',
    ]:
        rev_comp = biopsy.reverse_complement( seq )
        print seq, rev_comp
        if biopsy.reverse_complement( rev_comp ) != seq:
            raise RuntimeError, 'Double reverse complement failed'






def test_lcs():
    seqs = [
            # 'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagca',
            # 'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagca',
            'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagcagactgcatgcggacgggaccgcggatcgcgcggtgcgcctcagtgtacttccgaacgaatgagtcattaatagagcgctatatcgtaactgtctttgacgaagtataccgaaaccgtgcagccagacgtgatccgggcgttgtaaaggcgatcagcgccctaggagtaccatttttgccgtaggcttgcgtctcaaagaccagctggggcgtggtatcactcgtcagtacgatttctgccagatagatagcatagactgaaccttaggcccaatagggacacaattacccgagtgactgactggtctaaggggagtccccccttaaaacgttttacgtaatagcgggctccagaagcaaagcatcggtttgagccccagtactaaacgtttgagtgtttgctctcgtctgataggtaaaccgacaagagaaccaagctcaaggcgcggtaggtgcgccttgcgaactgttgatgccgtgagcgccaccatcccgtgcatcataggcagggagagaagaccacatggccttgcgaccgtatgagctgtttcagattaaatgccaacgggcatggtcggtgtccagcattttttgcagtcagctggtggtacacagtggggacaagaacgcctctggtagatgtcttctgaaggagtaactcatttcgttgaatcgaccttcccttgcgcttgaacgcggacctctagtctctctcgcagactggggtcgaaaatcaaggtagatatggaatgttccgcatgagggtagcgaccggatcgggcgtcaagtatatcctccctgctacgtccccctactagcctcagtccgcctcgaacctaggaagattggccacatcagcttggtggatgcctggtccatacttcagacccgagaatgttagacaggaccccatttggctcctttacgtacgatctatgtagacgcagtga',
            'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagcagactgcatgcggacgggaccgcggatcgcgcggtgcgcctcagtgtacttccgaacgaatgagtcattaatagagcgctatatcgtaactgtctttgacgaagtataccgaaaccgtgcagccagacgtgatccgggcgttgtaaaggcgatcagcgccctaggagtaccatttttgccgtaggcttgcgtctcaaagaccagctggggcgtggtatcactcgtcagtacgatttctgccagatagatagcatagactgaaccttaggcccaatagggacacaattacccgagtgactgactggtctaaggggagtccccccttaaaacgttttacgtaatagcgggctccagaagcaaagcatcggtttgagccccagtactaaacgtttgagtgtttgctctcgtctgataggtaaaccgacaagagaaccaagctcaaggcgcggtaggtgcgccttgcgaactgttgatgccgtgagcgccaccatcccgtgcatcataggcagggagagaagaccacatggccttgcgaccgtatgagctgtttcagattaaatgccaacgggcatggtcggtgtccagcattttttgcagtcagctggtggtacacagtggggacaagaacgcctctggtagatgtcttctgaaggagtaactcatttcgttgaatcgaccttcccttgcgcttgaacgcggacctctagtctctctcgcagactggggtcgaaaatcaaggtagatatggaatgttccgcatgagggtagcgaccggatcgggcgtcaagtatatcctccctgctacgtccccctactagcctcagtccgcctcgaacctaggaagattggccacatcagcttggtggatgcctggtccatacttcagacccgagaatgttagacaggaccccatttggctcctttacgtacgatctatgtagacgcagtga',
            # 'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagcagactgcatgcggacgggaccgcggatcgcgcggtgcgcctcagtgtacttccgaacgaatgagtcattaatagagcgctatatcgtaactgtctttgacgaagtataccgaaaccgtgcagccagacgtgatccgggcgttgtaaaggcgatcagcgccctaggagtaccatttttgccgtaggcttgcgtctcaaagaccagctggggcgtggtatcactcgtcagtacgatttctgccagatagatagcatagactgaaccttaggcccaatagggacacaattacccgagtgactgactggtctaaggggagtccccccttaaaacgttttacgtaatagcgggctccagaagcaaagcatcggtttgagccccagtactaaacgtttgagtgtttgctctcgtctgataggtaaaccgacaagagaaccaagctcaaggcgcggtaggtgcgccttgcgaactgttgatgccgtgagcgccaccatcccgtgcatcataggcagggagagaagaccacatggccttgcgaccgtatgagctgtttcagattaaatgccaacgggcatggtcggtgtccagcattttttgcagtcagctggtggtacacagtggggacaagaacgcctctggtagatgtcttctgaaggagtaactcatttcgttgaatcgaccttcccttgcgcttgaacgcggacctctagtctctctcgcagactggggtcgaaaatcaaggtagatatggaatgttccgcatgagggtagcgaccggatcgggcgtcaagtatatcctccctgctacgtccccctactagcctcagtccgcctcgaacctaggaagattggccacatcagcttggtggatgcctggtccatacttcagacccgagaatgttagacaggaccccatttggctcctttacgtacgatctatgtagacgcagtga',
            ]
    hits = biopsy.HitsVec()
    for s in seqs:
        hits.append( biopsy.analyse( s, 0.03 ) )
    lcs = biopsy.longest_common_subsequence( hits )
    print lcs



def test_pssm_score():
    # 'V$AP1_Q2'
    pssm_acc = biopsy.get_transfac_pssm_accession( 'V$DEAF1_01' );
    pssm_info = biopsy.get_pssm( pssm_acc )
    # print pssm_info.pssm
    seq = 'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagcagactgcatgcggacgggaccgcggatcgcgcggtgcgcctcagtgtacttccgaacgaatgagtcattaatagagcgctatatcgtaactgtctttgacgaagtataccgaaaccgtgcagccagacgtgatccgggcgttgtaaaggcgatcagcgccctaggagtaccatttttgccgtaggcttgcgtctcaaagaccagctggggcgtggtatcactcgtcagtacgatttctgccagatagatagcatagactgaaccttaggcccaatagggacacaattacccgagtgactgactggtctaaggggagtccccccttaaaacgttttacgtaatagcgggctccagaagcaaagcatcggtttgagccccagtactaaacgtttgagtgtttgctctcgtctgataggtaaaccgacaagagaaccaagctcaaggcgcggtaggtgcgccttgcgaactgttgatgccgtgagcgccaccatcccgtgcatcataggcagggagagaagaccacatggccttgcgaccgtatgagctgtttcagattaaatgccaacgggcatggtcggtgtccagcattttttgcagtcagctggtggtacacagtggggacaagaacgcctctggtagatgtcttctgaaggagtaactcatttcgttgaatcgaccttcccttgcgcttgaacgcggacctctagtctctctcgcagactggggtcgaaaatcaaggtagatatggaatgttccgcatgagggtagcgaccggatcgggcgtcaagtatatcctccctgctacgtccccctactagcctcagtccgcctcgaacctaggaagattggccacatcagcttggtggatgcctggtccatacttcagacccgagaatgttagacaggaccccatttggctcctttacgtacgatctatgtagacgcagtga'
    for i in range( len( seq ) - len( pssm_info.pssm )  + 1 ):
        s = biopsy.score_pssm( pssm_info.pssm, seq[i:] )
        p_binding = biopsy.get_p_binding(
                biopsy.get_odds_ratio(
                        s,
                        pssm_info.get_dist( True, False ),
                        pssm_info.get_dist( False, False ) ) )
        if p_binding > 0.05:
            print i, s, p_binding
    result = biopsy.HitVec()
    p_binding = biopsy.score_pssm_on_sequence( pssm_acc, seq, 0.05, result )
    print 'Got', len( result ), 'hits from', len( seq ), 'bases'
    print p_binding
# test_pssm_score()


def test_score_pssms():
    # 'V$AP1_Q2'
    transfac_pssms = biopsy.get_transfac_pssm_accessions( biopsy.get_default_transfac_pssm_filter() )
    print 'Got', len( transfac_pssms ), 'pssms'
    seq = 'tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagcagactgcatgcggacgggaccgcggatcgcgcggtgcgcctcagtgtacttccgaacgaatgagtcattaatagagcgctatatcgtaactgtctttgacgaagtataccgaaaccgtgcagccagacgtgatccgggcgttgtaaaggcgatcagcgccctaggagtaccatttttgccgtaggcttgcgtctcaaagaccagctggggcgtggtatcactcgtcagtacgatttctgccagatagatagcatagactgaaccttaggcccaatagggacacaattacccgagtgactgactggtctaaggggagtccccccttaaaacgttttacgtaatagcgggctccagaagcaaagcatcggtttgagccccagtactaaacgtttgagtgtttgctctcgtctgataggtaaaccgacaagagaaccaagctcaaggcgcggtaggtgcgccttgcgaactgttgatgccgtgagcgccaccatcccgtgcatcataggcagggagagaagaccacatggccttgcgaccgtatgagctgtttcagattaaatgccaacgggcatggtcggtgtccagcattttttgcagtcagctggtggtacacagtggggacaagaacgcctctggtagatgtcttctgaaggagtaactcatttcgttgaatcgaccttcccttgcgcttgaacgcggacctctagtctctctcgcagactggggtcgaaaatcaaggtagatatggaatgttccgcatgagggtagcgaccggatcgggcgtcaagtatatcctccctgctacgtccccctactagcctcagtccgcctcgaacctaggaagattggccacatcagcttggtggatgcctggtccatacttcagacccgagaatgttagacaggaccccatttggctcctttacgtacgatctatgtagacgcagtga'
    # seq = 'acatcat'
    # seq = 'gat'
    # hits = biopsy.HitVec()
    hits = biopsy.score_pssms_on_sequence(
            transfac_pssms,
            seq,
            0.05 )
    print hits
    print 'score_pssm_on_sequence: Got', len( hits ), 'hits from', len( seq ), 'bases'
    hits = biopsy.analyse(
            seq,
            0.05)
    # print hits
    print 'analyse: Got', len( hits ), 'hits from', len( seq ), 'bases'
# test_score_pssms()





def test_index():
    index = biopsy.MultiDimIdx( [ 4, 4, 4 ] )
    while index.increment():
        # print index.index, 'has zero index:', index.has_zero_index()
        pass






def test_lcs():
    lcs_test_seqs = [
                    [
                            'ACA',
                            'ACCA',
                    ],
                    [
                            'a1b2c3d4e',
                            'zz1yy2xx3ww4vv',
                    ],
                    [
                            'abcdgh',
                            'aedfhr',
                    ],
                    [
                            'abcdefghijklmnopqrstuvwxyz',
                            'a0b0c0d0e0f0g0h0i0j0k0l0m0n0o0p0q0r0s0t0u0v0w0x0y0z0',
                    ],
                    [
                            'abcdefghijklmnzyxwvutsrqpo',
                            'opqrstuvwxyzabcdefghijklmn',
                    ],
                    [
                            'TCTCAGTATGCGCACCTCTCCCGAGACACCTTTGAAGTC',
                            'GCCAGGTATGTGCATCTCCCCTCAATCTGCCTTTGAAGC',
                            'TATGGGGTATGTGCATCTTCCTCACCTGCCTTTGAAGCC',
                    ],
            ]

    for test_seq in lcs_test_seqs:
        cl = time.clock()
        lcs = \
                biopsy.LCS(
                        test_seq,
                        biopsy.get_char_for_nucleo,
                        biopsy.get_constant_score )
        print 'Took:', time.clock() - cl
        # print 'Strings:', lcs.strings
        print lcs.get_best()
        print

# profile.run( 'test_lcs()', 'lcs.prof' )
# sys.exit()


def test_analyse():
    sequences = {
            'mouse'         : 'TCTCAGTATGCGCACCTCTCCCGAGACACCTTTGAAGTCCTGGGTCCTTTGATGTGAGTGGGGAATGCGGTTTGCGGGGAGATCCTGAGACCCCTGTGGTCAGTGGCCCATAGAGCATGGGTGCTTGCTG',
            'dog'           : 'GCCAGGTATGTGCATCTCCCCTCAATCTGCCTTTGAAGCCCGGGGCCCTTTGATGTGAGTGGGGGATGCACCAGCCTGGGATCTTGCCAGCCTTGTGGTCAGTAGCCAACAGGGCACGGGGCTTGGTGGAGTTTGGAAATAAGAG',
            'human'         : 'TATGGGGTATGTGCATCTTCCTCACCTGCCTTTGAAGCCTGGGGCCCTTTGATGTTAGTGGGGAATGCACCAGCCTGGGGCTCTTGAGAGCTCTGCGGTCAGTGGCCAACAAGGAACAGGGACTGGGGAG',
    }

    print 'Analysing'
    hits = { }
    for name, seq in sequences.iteritems():
        hits[ name ] = biopsy.analyse( seq, 0.01 )

        print name, 'had', len( hits[ name ] ), 'hit(s)'
        print hits[ name ]



    for keys in [
                    [ 'mouse', 'dog' ],
                    [ 'mouse', 'human' ],
                    [ 'human', 'dog' ],
                    [ 'mouse', 'human', 'dog' ],
            ]:
        print keys
        seq, score = \
                biopsy.LCS(
                        [ hits[ key ] for key in keys ],
                        biopsy.get_char_for_hit,
                        biopsy.get_score_for_hit ).get_best()
        print score
        print ",".join( seq )
        print

if __name__ == '__main__':
    unittest.main()
