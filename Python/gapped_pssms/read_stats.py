#
# Copyright John Reid 2008,2009
#

"""
Reads ROC statistics.
"""


import sys, logging, pylab as P
from cookbook.dicts import dict_of
from infpy.roc import RocCalculator, plot_roc_points, label_plot, plot_random_classifier
from glob import glob

def load_stat_files(files, stats):
    for stats_file in files:
        logging.info('Loading statistics from %s', stats_file)
        f = open(stats_file)
        f.next()
        for line in f:
            bg, dataset, num_pssms, p_binding, tp, tn, fp, fn = line.split(';')
            if dataset not in datasets_to_ignore:
                datasets.add(dataset)
                bg_types.add(bg)
                num_pssms = int(num_pssms)
                p_binding = float(p_binding)
                tp = int(tp)
                tn = int(tn)
                fp = int(fp)
                fn = int(fn)
                stats[bg][num_pssms][dataset][p_binding] = RocCalculator(tp, fp, tn, fn)
                stats[bg][num_pssms]['overall'][p_binding] += RocCalculator(tp, fp, tn, fn)

logging.basicConfig(level=logging.INFO)
logging.getLogger('').addHandler(logging.FileHandler('read_stats.log'))
logging.info('Command line: %s', ' '.join(sys.argv))

# indexed by method, bg-type, num pssms, dataset, p(binding)
overall_stats = dict_of(dict_of(dict_of(dict_of(dict_of(RocCalculator)))))()
datasets = set(('overall',))
bg_types = set()
methods = []
datasets_to_ignore = set(('T99004',))
show_title = False

sys.argv.pop(0) # remove script name
while sys.argv:
    method = sys.argv.pop(0)
    methods.append(method)
    file_glob = sys.argv.pop(0)
    logging.info('Method: %s; glob: %s', method, file_glob)
    load_stat_files(glob(file_glob), overall_stats[method])

styles = [ '-', '--' ]
for bg in bg_types:
    for dataset in datasets:
        P.figure()
        for style, method in zip(styles, methods):
            stats = overall_stats[method]
            rocs = list(stats[bg][0][dataset].iteritems())
            rocs.sort()
            rocs = [r[1] for r in rocs]
            plot_roc_points(rocs, linestyle=style, color='k', label=method)
        label_plot()
        P.legend(loc='lower right')
        if show_title:
            P.title('%s - %s' % (dataset, bg))
        plot_random_classifier(label='random')
        P.savefig('ROC-%s-%s.png' % (bg, dataset))
        P.savefig('ROC-%s-%s.eps' % (bg, dataset))
        P.clf()
