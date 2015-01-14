#!/usr/bin/python

#   test-against-scipy.py for pairwise Version 1.0
#
#   A quick-and-dirty Python script to test if, for a given set of N random
#   points, the array of calculated pairwise distances output by pairwise is
#   the same as that output by scipy.
#
#   Note that here NumPy's allclose() is used to compare the arrays element-
#   wise, which returns True if compared elements have the same value to within
#   some very small error threshold. Some such error must be tolerated in these
#   comparisons because pairwise and scipy use different libraries to carry out
#   floating-point operations, leading to occasional small differences in the
#   results of otherwise identical calculations at high precisions.
#
#   Usage: python test.py
#
#   Returns 0 on success; returns 1 and prints error message on failure.

import sys
import os

sys.path.append(os.path.split(sys.path[0])[0])

from random import random
import numpy

try:
    import scipy.spatial
except Exception:
    exit(2)
    
import pairwise

n = 5000

p = [[random(), random(), random()] for i in xrange(n)]

r_s = scipy.spatial.distance.pdist(p)

r_p = pairwise.distances(p)

if not numpy.allclose(r_s, r_p):
    print("Error: outputs from pairwise and scipy are different.")
    exit(1)
