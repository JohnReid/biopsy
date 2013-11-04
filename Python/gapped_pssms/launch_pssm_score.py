#
# Copyright John Reid 2009
#

"""
Code to score the output of a algorithm using the test harness.
"""

import os, subprocess, logging, sys, Queue, threading, glam2
from optparse import OptionParser
from test_harness_2 import TestHarness

#
# Initialise the logging
#
format='%(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=format)
file_handler = logging.FileHandler('launch_pssm_score.log')
file_handler.setFormatter(logging.Formatter(format))
logging.getLogger('').addHandler(file_handler)
logging.info('Command line: %s', ' '.join(sys.argv))
logging.info('Current working directory: %s', os.getcwd())


#
# Parse the options
#
option_parser = OptionParser()
TestHarness.add_options(option_parser)
option_parser.add_option(
  "-j",
  "--num-threads",
  dest="num_threads",
  type='int',
  default=5,
  help="How many threads to run."
)
option_parser.add_option(
  "--pssm-dir",
  dest="pssm_dir",
  default='.',
  help="Directory PSSM files are stored in."
)
option_parser.add_option(
  "--glam2-format",
  dest="glam2_format",
  default=False,
  action="store_true",
  help="Algorithm output files are in GLAM2 format."
)
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%30s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)
if 1 != len(args):
    print >> sys.stderr, 'USAGE: %s <method>' % sys.argv[0]
    sys.exit(-1)
method = args[0]
logging.info('Method=%s', method)

harness = TestHarness(options)

def worker():
    while True:
        if options.glam2_format:
            pssm_file = q.get()
            fragment, fold, seed = glam2.interpret_output_filename(pssm_file)
        else:
            fragment, fold = q.get()
            pssm_file = os.path.join(options.pssm_dir, '%s-%d.pssm' % (fragment, fold))
        try:
            def execute(args):
                logging.info('Executing: %s', ' '.join(args))
                subprocess.check_call(args)

            logging.info('Scoring output: %s: %s, %d', pssm_file, fragment, fold)
            script_dir = os.path.dirname(__file__)
            args = [
                'python2.5',
                os.path.join(script_dir, 'test_score_pssm.py'),
                '-q',
                '-d', harness.options.data_dir, # where sequences are stored
                '-r', options.results_dir, # where to put output
                method,
                pssm_file,
                fragment,
                str(fold)
            ]
            if options.glam2_format:
                args.append('--glam2-format')

            # do positive dataset
            execute(args)

            # do negative datasets
            args.append(None)
            for bg in harness.options.backgrounds:
                args[-1] = bg
                execute(args)
        except:
            type, value, traceback = sys.exc_info()
            logging.exception('Exception (%s) caught: %s', type, value)
            sys.exc_clear()
        q.task_done()


#
# Create queue and worker threads
#
q = Queue.Queue()
for i in range(options.num_threads):
    t = threading.Thread(target=worker)
    t.setDaemon(True)
    t.start()


#
# Populate queue
#
if options.glam2_format:
    for best in glam2.best_files(options.pssm_dir):
        q.put(best)
else:
    for fragment in harness.options.fragments:
        for fold in harness.folds():
            q.put((fragment, fold))


#
# Wait until tasks are done
#
q.join()
