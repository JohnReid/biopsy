#!/usr/bin/env python2

import cPickle
import numpy
from rpy2.robjects import FloatVector, DataFrame
from rpy2.robjects.lib import ggplot2
from rpy2.robjects.packages import importr

pscores = cPickle.load(open('bootstrap.pickle'))
pscores.sort()
proportion = numpy.linspace(1, len(pscores), len(pscores)) / len(pscores)
dataf = DataFrame({
    'pscore': FloatVector(pscores),
    'proportion': FloatVector(proportion),
})

grdevices = importr('grDevices')
#grdevices.postscript(file="pscores.eps", width=512, height=512)
grdevices.postscript(file='pscores.eps')
(
    ggplot2.ggplot(dataf)
    + ggplot2.aes_string(y='pscore', x='proportion')
    + ggplot2.geom_point()
    + ggplot2.scale_x_log10()
    + ggplot2.scale_y_log10()
    + ggplot2.stat_smooth(method='lm')
).plot()
grdevices.dev_off()
