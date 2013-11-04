#
# Copyright John Reid 2009
#

"""
Code to launch the GLAM2 algorithm on the cross-validation of several fragments concurrently.
"""

import os, subprocess, logging, sys, Queue, threading
from optparse import OptionParser
from glam2 import *

#
# Initialise the logging
#
format='%(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=format)
file_handler = logging.FileHandler('launch_glam2_cross_validate.log')
file_handler.setFormatter(logging.Formatter(format))
logging.getLogger('').addHandler(file_handler)
logging.info('Command line: %s', ' '.join(sys.argv))
logging.info('Current working directory: %s', os.getcwd())


#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-j",
  "--num-threads",
  dest="num_threads",
  type='int',
  default=5,
  help="How many threads to run."
)
option_parser.add_option(
  "-n",
  dest="n",
  type='int',
  default=5000000,
  help="How many iterations to wait without improvement before stopping."
)
option_parser.add_option(
  "-f",
  "--num-folds",
  dest="num_folds",
  type='int',
  default=5,
  help="How many folds in the cross validation."
)
option_parser.add_option(
  "-s",
  "--num-seeds",
  dest="num_seeds",
  type='int',
  default=8,
  help="How many seeds to use."
)
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%30s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)

data_dir = args.pop(0)
fragments = args
logging.info('Data directory: %s' % data_dir)
logging.info('Fragments: %s' % ' '.join(fragments))

def worker():
    while True:
        fragment, cross_fold_index, seed = q.get()
        try:
            tag = make_tag(fragment, cross_fold_index, seed)
            fasta = os.path.join(data_dir, '%strimRM-train-x%d.fa' % (fragment, cross_fold_index))
            args = [
                'glam2', '-2', '-r', '1', '-I', '.1', '-J', '99999.0',
                '-n', str(options.n),
                '-s', str(seed),
                'n', fasta
            ]
            logging.info('%s: Executing: %s', tag, ' '.join(args))
            output_file = output_filename(fragment, cross_fold_index, seed)
            logging.info('%s: Sending output to: %s', tag, output_file)
            output = open(output_file, 'w')
            process = subprocess.Popen(args, stdout=output)
            process.wait()
            output.close()
            retcode = process.returncode
            if 0 != retcode:
                logging.error('%s: Exit status: %d: %s', tag, retcode, args)
        except:
            type, value, traceback = sys.exc_info()
            logging.error('%s: Exception (%s) caught: %s', tag, type, value)
            sys.exc_clear()
        q.task_done()


q = Queue.Queue()
for i in range(options.num_threads):
    t = threading.Thread(target=worker)
    t.setDaemon(True)
    t.start()

for fragment in fragments:
    for cross_fold_index in xrange(1, options.num_folds+1):
        for seed in xrange(1, options.num_seeds+1):
            q.put((fragment, cross_fold_index, seed))


q.join()       # block until all tasks are done
