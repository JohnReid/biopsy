
import string_kernel, numpy

def mismatch_weight( a, b ):
    import math
    assert isinstance( a, str )
    assert isinstance( b, str )
    if a == b:
        return math.log( 1.0 )
    else:
        return math.log( 0.5 )

def test_string_kernel():
    mutations = numpy.eye( 4, dtype = numpy.float64 )
    c = string_kernel.MismatchKernelCalculator.create( mutations )
    seqs = string_kernel.convert_to_sequence_vec(
            [
                    [ 1, 1, 2, 3 ],
                    [ 0, 1, 2, 3 ],
                    [ 0, 1, 2, 3 ],
                    [ 2, 1, 2, 2 ],
            ]
    )
    k = c.calculate( seqs )
    # print k

    strings = [
            'a',
            't',
            'a',
            ''
    ]
    k, sym_to_idx, idx_to_sym = string_kernel.mismatch_kernel( strings, mismatch_weight )
    print k
    print sym_to_idx
    print idx_to_sym


import Bio.SubsMat.MatrixInfo
blosum62 = Bio.SubsMat.MatrixInfo.blosum62
sub_mat = Bio.SubsMat.SeqMat(
        data = blosum62,
        mat_type = Bio.SubsMat.LO,
        mat_name = 'BLOSUM62'
)

class SubMatWeight( object ):
    def __init__( self, sub_mat ):
        self.sub_mat = sub_mat
    def __call__( self, a, b ):
        k = a < b and (a,b) or (b,a)
        if not self.sub_mat.has_key( k ):
            raise RuntimeError(
                    'Substitution matrix does not have key: %s' % str( k )
            )
        return self.sub_mat[ k ]

k, sym_to_idx, idx_to_sym = string_kernel.mismatch_kernel(
        [
                'ABCDEFGHIKLMNPQRSTVWXYZ',
                'AAAE',
                'ABCDEFGHIKLMNPQRSTVWXYZ',
        ],
        SubMatWeight( sub_mat )
)

def load_scop(
        scop_data_dir = 'c:/data/scop',
        scop_version = '1.71'
):
    """Return the entire scop hierachy parsed from the data files"""
    import corebio.resource.scop
    return corebio.resource.scop.Scop.parse_files(
            cla_file = open( '%s/dir.cla.scop.txt_%s' % (scop_data_dir, scop_version), 'r' ),
            des_file = open( '%s/dir.des.scop.txt_%s' % (scop_data_dir, scop_version), 'r' ),
            hie_file = open( '%s/dir.hie.scop.txt_%s' % (scop_data_dir, scop_version), 'r' )
    )
if None == scop: scop = load_scop()
