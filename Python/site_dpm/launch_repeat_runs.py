#
# Copyright John Reid 2009
#


"""
Code to launch the site DPM framework several times concurrently.
"""


import os, subprocess, logging, sys, Queue, threading, time, random


def ensure_dir_exists(dir):
    "Makes a directory if it does not already exist."
    if not os.access(dir, os.X_OK):
        logging.info('Making directory: %s', dir)
        os.makedirs(dir)


def get_new_dir(template_dir, run_index):
    return os.path.join(template_dir, '%04d' % run_index)


def make_link(src, dst):
    if os.path.exists(dst):
        os.unlink(dst)
    os.symlink(src, dst)


def make_new_dir_from_template(template_dir, new_dir):
    "Copies or links data from template directory into new directory."
    # make sure directories exist
    ensure_dir_exists(os.path.join(new_dir, 'cache'))

    for ext in [
        'hit-counts',
        'ucsc-analysis',
    ]:
        # link cache file
        cached_filename = '%s.%s' % (tag, ext)
        if os.path.exists(os.path.join(template_dir, 'cache', cached_filename)):
            template_filename = os.path.join('..', '..', 'cache', cached_filename)
            new_filename = os.path.join(new_dir, 'cache', cached_filename)
            make_link(template_filename, new_filename)

    # link options file if it exists
    if os.path.exists(os.path.join(template_dir, 'options.py')):
        make_link(os.path.join('..', 'options.py'), os.path.join(new_dir, 'options.py'))


def execute_cmd(args, cwd):
    "Execute the command."
    logging.info('Executing: %s in %s', args, cwd)
    retcode = subprocess.call(args, cwd=cwd)
    if 0 != retcode:
        logging.error('Exit status: %d: %s', retcode, args)


def worker():
    "Worker thread to execute tasks from the queue."
    while True:
        run_index = q.get()
        logging.info('Run #: % 3d' % run_index)
        new_dir = get_new_dir(template_dir, run_index)
        make_new_dir_from_template(template_dir, new_dir)
        args = [
            'python',
            '../../../tp_go_analysis.py',
            tag
        ]
        execute_cmd(args, new_dir)
        args = [
            'python',
            '../../../validate.py',
            tag
        ]
        execute_cmd(args, new_dir)
        q.task_done()


#
# Set up logging
#
logging.basicConfig(level=logging.DEBUG)
file_handler = logging.FileHandler('launch-repeat-runs.log')
file_handler.setFormatter(logging.Formatter('%(asctime)s:%(levelname)s:%(message)s'))
logging.getLogger().addHandler(file_handler)


#
# Get positional arguments
#
try:
    logging.info('Command line: %s', ' '.join(sys.argv))
    script_name = sys.argv.pop(0)

    template_dir = sys.argv.pop(0)
    logging.info('Template directory: %s' % template_dir)

    num_runs = int(sys.argv.pop(0))
    logging.info('# runs: %d' % num_runs)

    tag = sys.argv.pop(0)
    logging.info('Tag: %s' % tag)

except:
    logging.error('USAGE: python %s <template directory> <number of runs> <tag>', script_name)
    raise


#
# Set up worker threads
#
num_worker_threads = 4
stagger = 3
logging.info('Using %d worker threads' % num_worker_threads)

logging.info('Initialising worker threads')
q = Queue.Queue()
for i in range(num_worker_threads):
    t = threading.Thread(target=worker)
    t.setDaemon(True)
    t.start()

logging.info('Populating task queue')
for run_index in xrange(num_runs):
    q.put(run_index)
    time.sleep(1.)


#
# Wait until completed
#
logging.info('Blocking until all tasks completed')
q.join()       # block until all tasks are done
logging.info('All tasks completed')
