#!/usr/bin/python

#   time_scipy.py for pairwise Version 1.0
#
#   A quick-and-dirty Python script which returns the time taken for scipy to
#   calculate all pairwise distances for a set of POINTS randomly generated
#   points. This is provided as a point of reference against which to benchmark
#   pairwise.
#
#   Usage: python time-scipy.py POINTS
#
#       with POINTS a positive integer.

import sys
import timeit

n = int(sys.argv[1])

setup = list()
setup.append("import scipy.spatial")
setup.append("from random import random")
setup.append("points = [[random(), random(), random()] for i in xrange(%d)]" % n)
setup = "; ".join(setup)

t = timeit.Timer("scipy.spatial.distance.pdist(points)", setup)
print(t.timeit(1))
