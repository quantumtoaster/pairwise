#!/usr/bin/python

#   time_pairwise.py for pairwise Version 1.0
#
#   A quick-and-dirty Python script which returns the time taken for pairwise
#   to calculate all pairwise distances for a set of POINTS randomly generated
#   points. pairwise will use THREADS threads; THREADS = 0 means try to detect
#   the number of processors made available by the host and then use that many
#   threads, or use one thread if the number of processors cannot be found. The
#   verbosity level is set by VERBOSE; VERBOSE = 0 means don't print anything
#   to stdout; VERBOSE = 1 means print some information about each worker
#   thread as it is launched.
#
#   Usage: python time-pairwise.py POINTS [THREADS = 1] [VERBOSE = 0]
#
#       with POINTS, THREADS and VERBOSE all positive integers.

import sys
import timeit

n = int(sys.argv[1])

if len(sys.argv) > 2:
    t = int(sys.argv[2])
else:
    t = 1
    
if len(sys.argv) > 3:
    v = int(sys.argv[3])
else:
    v = 0

setup = list()
setup.append("import sys")
setup.append("import os")
setup.append("sys.path.append(os.path.split(sys.path[0])[0])")
setup.append("import pairwise")
setup.append("from random import random")
setup.append("p = [[random(), random(), random()] for i in xrange(%d)]" % n)
setup = "; ".join(setup)

t = timeit.Timer("pairwise.distances(points = p, threads = %d, verbose = %d)" % (t, v), setup)
print(t.timeit(1))
