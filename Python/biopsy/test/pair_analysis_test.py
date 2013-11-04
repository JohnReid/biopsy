#
# Copyright John Reid 2006
#

import biopsy, unittest

class TestMultinomialSample(unittest.TestCase):
    def test_samples( self ):
        for y in [
                [10]*3,
                [10]*4,
                [28,1,1],
                [28,1,1,1],
                [1,1,28],
                [0,0,30],
                [16,7,7],
                [20,5,5,5],
                [16,4,4,4],
                [12,3,3,3],
                [8,2,2,2],
                [4,1,1,1],
                [4,1,1,1,1],
        ]:
            print y
            m = biopsy.MultinomialSample( y )
            ll_uniform = m.ll_under_multinomial()
            ll_dirichlet = m.ll_under_dirichlet()
            lor = ll_uniform - ll_dirichlet
            model_cmp = m.model_comparison()
            print ll_uniform
            print ll_dirichlet
            print 'Log odds ratio for uniform: %f - %s' % ( ll_uniform - ll_dirichlet, y )
            print 'Lor      : %f' % ( lor )
            print 'Model cmp: %f' % ( model_cmp )

if __name__ == '__main__':
    unittest.main()
