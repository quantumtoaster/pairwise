#!/usr/bin/python

#   profile.py for pairwise Version 1.0
#
#   A quick-and-dirty Python script to profile the performance of pairwise both
#   in single- and multi-threaded mode against that of scipy. The time taken to
#   calculate all pairwise distances for a set of N randomly generated points
#   is recorded REPS times using pairwise in multi-threaded mode (using MULTI
#   threads), using pairwise in single-threaded mode, and then using scipy. The
#   mean average time taken for each such run is then appended to a date- and
#   time-stamped report file in the profiles directory.
#
#   Usage: python profile.py
#
#   REPS and MULTI below can be adjusted as desired.

import sys
import os
import datetime
import subprocess

if __name__ == "__main__":

    sys.path.append(os.path.split(sys.path[0])[0])

    try:
        os.mkdir("profiles")
    except:
        pass

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    out = open("profiles/%s_profile.csv" % timestamp, "wb")
    
    REPS = 100
    
    MULTI = 24

    runs = [10, 50, 100, 500, 1000, 3792, 5000, 10000, 20000, 30000, 40000, 50000]
    
    for n in runs:
        
        times1 = list()
        times2 = list()
        times3 = list()
    
        for i in xrange(REPS):
            
            o1 = subprocess.Popen(["%s %d %d 0" % (os.path.join(sys.path[0], "time-pairwise.py"), n, MULTI)], stdout = subprocess.PIPE, shell = True)
            times1.append(float(o1.communicate()[0].strip()))
            
            o2 = subprocess.Popen(["%s %d %d 0" % (os.path.join(sys.path[0], "time-pairwise.py"), n, 1)], stdout = subprocess.PIPE, shell = True)
            times2.append(float(o2.communicate()[0].strip()))
            
            o3 = subprocess.Popen(["%s %d" % (os.path.join(sys.path[0], "time-scipy.py"), n)], stdout = subprocess.PIPE, shell = True)
            times3.append(float(o3.communicate()[0].strip()))
        
        out.write("%07d,pairwise,%02d,%e,%e,%e\n" % (n, MULTI, sum(times1) / REPS, min(times1), max(times1)))
        out.write("%07d,pairwise,01,%e,%e,%e\n" % (n, sum(times2) / REPS, min(times2), max(times2)))
        out.write("%07d,   scipy,01,%e,%e,%e\n" % (n, sum(times3) / REPS, min(times3), max(times3)))
        
        out.flush()
        os.fsync(out)
        
    out.close()

