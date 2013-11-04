#
# Copyright John Reid 2006
#

from _biopsy import *

class MultiDimIdx:
    """An index into a multi dimensional array
    """

    def __init__( self, dimensions, index = None ):
        self.dimensions = dimensions
        if None == index:
            self.index = [ 0 ] * len( self.dimensions )
        else:
            self.index = [ i for i in index ]

    def __repr__( self ):
        return ",".join( str(i) for i in self.index )

    def increment( self ):
        """Move the index along one. Returns false when hits end"""
        i = len( self.index ) - 1
        while -1 != i:
            self.index[ i ] += 1
            if self.dimensions[ i ] == self.index[ i ]:
                self.index[ i ] = 0
                i -= 1
            else:
                return True
        return False


    def has_zero_dimension( self ):
        """Is one of the dimensions empty?"""
        return 0 in self.dimensions

    def is_at_end( self ):
        """Did we reach the end?"""
        for i, j in enumerate( self.index ):
            if self.dimensions[ i ] <= j:
                return True
        return False


    def get_num_data( self ):
        """Number of entries in multi dim array"""
        result = 1
        for length in self.dimensions:
            result *= length
        return result


    def get_data_index( self ):
        """Returns single dimension index"""
        result = 0
        base = 1
        for i, j in enumerate( self.index ):
            result += base * j
            base *= self.dimensions[ i ]
        return result

    def get_last( self ):
        result = MultiDimIdx( self.dimensions )
        result.index = [ d - 1 for d in self.dimensions ]
        return result

    def has_zero_index( self ):
        for i in self.index:
            if 0 == i:
                return True
        return False

    def copy( self ):
        """Returns a copy of this index"""
        return MultiDimIdx( self.dimensions, self.index )

    def get_all_previous( self ):
        """Gets all indices one before this one in each dimension"""
        result = []
        for i in range( len( self.index ) ):
            if 0 != self.index[ i ]:
                index_copy = self.copy()
                index_copy.index[ i ] -= 1
                result.append( index_copy )
        return result

    def get_previous( self ):
        """Gets the index which is one smaller in every dimension"""
        return MultiDimIdx( self.dimensions, [ i - 1 for i in self.index ] )



















class LCS:
    """Implements a longest common subsequence algorithm

    Uses a simple dynamic programming approach
    """

    def get_score( self, index ):
        return self.table[ index.get_data_index() ][ 1 ]

    def get_subsequence( self, index ):
        return self.table[ index.get_data_index() ][ 0 ]

    def get_best_score_from( self, indices ):
        best_score = None
        best_index = None
        for prev in indices:
            prev_idx = prev.get_data_index()
            # print 'Prev:', prev.index
            prev_score = self.table[ prev_idx ][ 1 ]
            if None == best_score or best_score < prev_score:
                best_index = prev
                best_score = prev_score
        return ( best_score, best_index )

    def are_all_chars_same( self, index ):
        # are the characters all the same?
        c = None
        for i, j in enumerate( index.index ):
            if None == c:
                # print index
                c = self.get_char_fn( self.strings[ i ][ j ] )
            else:
                if c != self.get_char_fn( self.strings[ i ][ j ] ):
                    c = None
                    break
        return c

    def get_best( self ):
        index = MultiDimIdx( self.lengths ).get_last()
        if index.has_zero_dimension():
            return ( [ ], 0 )
        return self.table[ index.get_data_index() ]

    def __init__( self, strings, get_char_fn, get_score_fn ):

        self.get_char_fn = get_char_fn
        self.get_score_fn = get_score_fn

        # get the set containing characters in every string
        self.char_universe = None
        for string in strings:
            s = set( get_char_fn( element ) for element in string )
            if None == self.char_universe:
                self.char_universe = s
            else:
                self.char_universe = self.char_universe & s
        # print 'Universe is:', self.char_universe

        self.strings = []
        for string in strings:
            self.strings.append(
                    filter(
                            lambda c: self.get_char_fn( c ) in self.char_universe,
                            string ) )

        # the lengths of the strings
        self.lengths = [ len( s ) for s in self.strings ]
        # print 'Lengths are:', self.lengths
        print 'Running time proportional to ', reduce( lambda x, y: x * y, self.lengths )


        # an index into our table
        index = MultiDimIdx( self.lengths )

        # dynamic programming multi-dimensional array
        # each element is tuple of the best common subsequence up to that point
        # and its score
        self.table = [ None ] * index.get_num_data()

        # iterate over indices
        while not index.is_at_end():
            data_index = index.get_data_index()

            # find the best from previous indices, i.e. those that we can just carry over
            best_score, best_index = self.get_best_score_from( index.get_all_previous() )
            if None != best_index:
                best_seq = self.get_subsequence( best_index )
            else:
                best_seq = []

            # check that we don't have better now
            # are all the characters the same?
            c = self.are_all_chars_same( index )
            if None != c:
                # do we have a better subsequence?
                score = None
                # find the lowest score in all the strings
                for i, j in enumerate( index.index ):
                    s = self.get_score_fn( self.strings[ i ][ j ] )
                    if None == score or s < score:
                        score = s
                seq = []
                if not index.has_zero_index():
                    diagonal_index = index.get_previous()
                    diagonal_score = self.get_score( diagonal_index )
                    if None != diagonal_score:
                        score += diagonal_score
                    seq = [ x for x in self.get_subsequence( diagonal_index ) ]
                seq.append( c )
                if None == best_score or best_score < score:
                    best_seq = seq
                    best_score = score

            self.table[ data_index ] = ( best_seq, best_score )
            # print self.table[ MultiDimIdx( self.lengths, [ 1, 1 ] ).get_data_index() ]
            # print index.index, best_score, best_seq, c

            if not index.increment(): break
