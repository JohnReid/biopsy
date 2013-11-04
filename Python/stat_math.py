import scipy.special, math

def dirichlet_log_pdf( parameters, values ):
    result = sum = 0.0
    assert len( values ) == len( parameters )
    for v, p in zip( values, parameters ):
        print "%f,%f" % (p,v)
        result += (p-1) * math.log( v ) - scipy.special.gammaln( p )
        sum += p
    return result + scipy.special.gammaln( sum )

print dirichlet_log_pdf(
        [ 100 ] * 4,
        [ .25 ] * 4
)

print dirichlet_log_pdf(
        [ 100 ] * 4,
        [
                0.23979402951046996,
                0.27145845441922739,
                0.25459264428166417,
                0.23415487178863836
        ]
)
