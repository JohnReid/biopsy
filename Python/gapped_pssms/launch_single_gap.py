
#
# Copyright John Reid 2009
#

"""
Code to launch the single gap algorithm on several fragments concurrently.
"""

import os, subprocess, logging, sys, Queue, threading
from optparse import OptionParser

def get_code_dir():
    logging.info('Locating code directory.')
    candidate_dirs = [
        '/home/reid/Dev/MyProjects/Bio/Python/gapped_pssms',
        '/home/john/Dev/MyProjects/Bio/Python/gapped_pssms',
    ]
    for candidate in candidate_dirs:
        if os.path.exists(candidate):
            return candidate
    raise ValueError('Cannot find code directory.')

logging.basicConfig(level=logging.DEBUG)
logging.getLogger('').addHandler(logging.FileHandler('launch_single_gap.log'))
logging.info('Command line: %s', ' '.join(sys.argv))
sys.argv.pop(0)

fragment_pattern = sys.argv.pop(0)
logging.info('Fragment pattern: %s' % fragment_pattern)

output_dir = sys.argv.pop(0)
logging.info('Output directory: %s' % output_dir)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

num_worker_threads = int(sys.argv.pop(0))
logging.info('# subprocesses: %d' % num_worker_threads)

fragments = sys.argv.pop(0).split()
logging.info('Fragments: %s' % ' '.join(fragments))

single_gap_args = sys.argv
logging.info('Single gap args: %s' % ' '.join(single_gap_args))

code_dir = get_code_dir()

def worker():
    while True:
        fragment = q.get()
        logging.info('Analysing %s' % fragment)
        args = [
            'python',
            os.path.join(code_dir, 'run_single_gap.py'),
            '--fasta', fragment_pattern % fragment,
            '--tag', fragment,
            '--output', output_dir
        ]
        args += single_gap_args
        logging.info('Executing: %s', args)
        retcode = subprocess.call(args)
        if 0 != retcode:
            logging.error('Exit status: %d: %s', retcode, args)
        q.task_done()

q = Queue.Queue()
for i in range(num_worker_threads):
    t = threading.Thread(target=worker)
    t.setDaemon(True)
    t.start()

for item in fragments:
    q.put(item)

q.join()       # block until all tasks are done
