#
# Copyright John Reid 2008
#

import os, time

def timestamp():
    return time.strftime('%Y-%m-%d--%H-%M-%S')

def make_output_dir(base_dir):
    dir = os.path.join(base_dir, timestamp())
    os.makedirs(dir)
    return dir
