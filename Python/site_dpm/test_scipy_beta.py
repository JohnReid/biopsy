import scipy.stats as S, numpy as N
from rpy2.robjects import r

alpha, beta = 0.0710, 0.4222
for i in xrange(20):
    x_from_scipy = S.beta.rvs(alpha, beta)
    x_from_R = r.rbeta(1, alpha, beta)
    print 'Alpha=%.2e; Beta=%.2e; scipy.stats.beta.rvs=%.2e; R.rbeta=%.2e' % (alpha, beta, x_from_scipy, x_from_R[0])
    alpha /= 10.
    beta /= 10.
