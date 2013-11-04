#
# Copyright John Reid 2006
#

import biopsy.gapped_pssms, unittest, numpy

class TestSite( unittest.TestCase ):
    def test_site( self ):
        "test some random choices for Site"
        from numpy.random import randint, uniform
        for test_site in xrange( 100 ):
            start = randint( 100 )
            rev_comp = uniform() > 0.5
            length = randint( 3, 20 )
            site = biopsy.gapped_pssms.Site( start, rev_comp, length )
            for test in xrange( 20 ):
                model_idx = randint( site.length )
                assert( model_idx == site.to_model_idx( site.to_seq_coord( model_idx ) ) )
                seq_coord = randint( 100 )
                assert( seq_coord == site.to_seq_coord( site.to_model_idx( seq_coord ) ) )

    def test_gapped_site( self ):
        from numpy.random import randint, uniform
        for i in xrange( 100 ):
            length = 4
            for gap_idx in xrange( length - 1 ):
                site = biopsy.gapped_pssms.GappedSite(
                        start = 10,
                        rev_comp = False,
                        length = length,
                        gap_idx = gap_idx
                )
                for r in xrange( site.length ):
                    # print gap_idx, r, site.idx( r ), site.r( site.idx( r ) )
                    assert( r == site.r( site.idx( r ) ) )
                for idx in xrange( site.length ):
                    # print gap_idx, idx, site.r( idx ), site.idx( site.r( idx ) )
                    assert( idx == site.idx( site.r( idx ) ) )

    def test_r_and_offset( self ):
        K = 5
        for gap_position in xrange( K ):
            for has_gap in ( False, True ):
                for rev_comp in ( False, True ):

                    # test offset -> r -> offset
                    for offset in xrange( K + 1 ):
                        r = biopsy.gapped_pssms.r( K, gap_position, offset, has_gap, rev_comp )
                        calc_offset = biopsy.gapped_pssms.offset( K, gap_position, r, has_gap, rev_comp )
                        if calc_offset != calc_offset:
                            print ' K, gap position, has gap, rev comp, offset, r, calc offset'
                            error = '%2d,%13d,%8s,%9s,%7d,%2d,%12d' % (
                                    K,
                                    gap_position,
                                    str(has_gap),
                                    str(rev_comp),
                                    offset,
                                    r,
                                    calc_offset
                            )
                            print error
                            raise RuntimeError( error )

                    # test r -> offset -> r
                    for r in xrange( K + 1 ):
                        offset = biopsy.gapped_pssms.offset( K, gap_position, r, has_gap, rev_comp )
                        r = biopsy.gapped_pssms.r( K, gap_position, offset, has_gap, rev_comp )
                        calc_r = biopsy.gapped_pssms.r( K, gap_position, offset, has_gap, rev_comp )
                        if r != calc_r:
                            print ' K, gap position, has gap, rev comp, r, offset, calc r'
                            error = '%2d,%13d,%8s,%9s,%2d,%7d,%7d' % (
                                    K,
                                    gap_position,
                                    str(has_gap),
                                    str(rev_comp),
                                    r,
                                    offset,
                                    calc_r
                            )
                            print error
                            raise RuntimeError( error )


if __name__ == '__main__':
    unittest.main()
