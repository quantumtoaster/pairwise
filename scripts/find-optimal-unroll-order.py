#!/usr/bin/python

#   find-optimal-unroll-order.py for pairwise Version 1.0
#
#   A quick-and-dirty Python script to find the order of manual loop unrolling
#   in pairwise's distance()'s time-intensive loop, such that distance()
#   returns as quickly as possible on the host.
#
#   For time-critical applications of pairwise, run this script on the target
#   host - it will probably take a while. Then run build.py with the value of
#   --unroll set to the optimal unroll order reported.
#
#   Usage: python test.py
#
#   Returns 0 on success; returns 1 and prints error message on failure.

import os
import sys
import subprocess

if __name__ == "__main__":

    n = 5000
    t = 1
    
    reps = 1000
    
    sys.path.append(os.path.split(sys.path[0])[0])
    
    averages = list()
    
    sys.stdout.write("find-optimal-unroll-order.py: finding optimal unroll order.\n")
    sys.stdout.flush()
    
    i = 0
    
    while True:
    
        ret = os.system("python setup.py build_ext --inplace --unroll-order=%d" % i)
        
        if ret: exit(1)
        
        times = list()
        
        for j in xrange(reps):
        
            proc = subprocess.Popen(["%s %d %d 0" % (os.path.join(sys.path[0], "time-pairwise.py"), n, t)], stdout = subprocess.PIPE, shell = True)
            times.append(float(proc.communicate()[0].strip()))
        
        average = sum(times) / reps
        
        averages.append(average)
        
        sys.stdout.write("find-optimal-unroll-order.py: order %d, averaged %e s.\n" % (i, average))
        sys.stdout.flush()
        
        i += 1
        
        if len(averages) > 2:
            if averages[-1] > averages[-2] and averages[-2] > averages[-3]:
                break
        
    optimal = averages.index(min(averages))
    
    sys.stdout.write("find-optimal-unroll-order.py: optimal unroll order is %d.\n" % optimal)
    sys.stdout.flush()
