
    pairwise Version 1.0
    
    A Python extension module for the multi-threaded calculation of pairwise
    distances from a set of points in three-dimensional space.
    
    Copyright (C) 2015 Adam J. Marsden


    Overview
    ========
    
      pairwise is an extension module for Python 2.7 written in GNU C. It was
    developed and tested on hosts running Ubuntu 14.04 and Linux Mint 17, where
    it requires the python-dev and python-numpy packages to be built.
    
      pairwise was developed using NumPy version 1.8.2; compatibility with
    other versions is unknown.
    
      pairwise carries out the minimum number of calculations required to find
    the Euclidean distances between all pairs of input points.
    
      In my speed tests, running in single-threaded mode, pairwise's
    distances() performs slightly better than the SciPy equivalent
    (scipy.spatial.distance.pdist()) for large sets of input points (of greater
    than about 200) and up to an order of magnitude better for small sets (of
    less than about 200). Running in multi-threaded mode, pairwise's
    distances() can perform significantly better than the SciPy equivalent for
    sets of input points of greater than about 1000, and an order of magnitude
    better for sets greater than about 5000. Naturally comparative performance
    depends strongly on the host platform and, in multi-threaded mode, the
    chosen ratio of threads-to-processors.
    
      pairwise was written for my own use in structural biology, and also as a
    learning exercise in writing modules for Python in C. To this end, I
    welcome any advice on how I might make it better. I release pairwise under
    the GNU General Public License Version 2 in hopes that it might be of use
    to others.


    Getting Started Quickly
    =======================
    
        $ python setup.py build
        $ sudo python setup.py install
    
   Then in Python,
    
        import pairwise
        from random import random
        
        p = [[random(), random(), random()] for i in xrange(10000)]
        
        d = pairwise.distances(points = p, threads = 10)
        
    and to get the distance between points p[100] and p[-1],
    
        i = pairwise.index(10000, 100, -1)
        
        my_distance = d[i]


    Installation
    ============
    
      pairwise is fully contained within a single file, source/pairwise.c. A
    custom Python distutils script setup.py handles compilation and
    installation on the host in the usual way. Note that for compilation to
    succeed the host must make available the Python 2.7 and NumPy development
    headers. On hosts running Ubuntu and Linux Mint (and probably other Debian
    distributions too), first run,
    
        $ apt-get install python-dev python-numpy
    
    to get these if you don't have them already.
    
      pairwise can be built in-place using,
        
        $ python setup.py build_ext --inplace [--unroll-order=INT]
    
    where --unroll-order is non-standard and optional.
    
      setup.py manually unrolls the time-intensive calculation loop at the
    centre of pairwise's distances() a number of times as specified by
    --unroll-order. On some hosts I've seen some small speed increases by
    choosing a value for --unroll-order which is greater than zero. If you
    don't care about micro-optimisations feel free to leave --unroll-order
    unspecified, in which case no loop unrolling will be done. If you do choose
    to specify --unroll-order with a value which is greater than zero, setup.py
    will write an appropriately unrolled copy of the pairwise's source to
    source/pairwise-unrolled.c, and then that will be compiled.
    
      If you choose to build pairwise in-place and all goes well you should end
    up with a shared object file pairwise.so in the distribution root. To use
    pairwise in Python scripts without permanently installing it on the host
    simply copy pairwise.so into the same directory as the Python script which
    intends to use it and,
    
        import pairwise
    
    as usual. Alternatively add the location of pairwise.so to your PYTHONPATH
    environment variable.
    
      If, on the other hand, you want to install pairwise on the host
    permanently (such that it will be available to all Python scripts) run,
    
        $ python setup.py build [--unroll-order=INT]
        $ sudo python setup.py install
    
    where the behaviour of --unroll-order is the same as above.


    Methods
    =======
    
      pairwise provides two methods: distances() and index().
        
      distances() accepts a two-dimensional array of real signed numbers
    representing the positions of a set of points in three-dimensional space.
    This array must have form,
    
        points = [[x1, y1, z1], [x2, y2, z2], ...]
    
    where xi, yi and zi are the x-, y-, and z-coordinates of the ith point.
    Provided it is of that form, this array can be any set of appropriately
    nested Python objects, so long as its first and second dimensions support
    Python's sequence protocol. In practice this means that this array could be
    a list of lists of floats, a tuple of tuples of integers, or even a two-
    dimensional NumPy array of floats. Internally all numbers are converted to
    C doubles and, at this time, buffer overflow checking is not done.
    
      distances() further accepts two optional arguments. These are threads and
    verbose, and are both positive integers.
    
      threads is the number of threads over which to distribute the pairwise
    calculations to be done. Internal calculations ensure that this
    distribution is as even as possible. A special value of 0 for threads asks
    distances() to probe /proc/cpuinfo if it is made available by the host, and
    to count the number of "processor" entries therein. If this count is
    greater than zero, distances() will then use that many threads for its
    calculations. If, on the other hand, /proc/cpuinfo cannot be read or
    contains no "processor" entries, distances() will raise a Python exception.
    The default value for threads is 1.
    
      verbose is the level of verbosity with which to run distances(). Its
    truthiness is evaluated, so meaningful values are 0 or greater-than-one. If
    verbose is set to 0, distances() executes silently. Otherwise if verbose is
    set to greater-than-one, distances() prints a small number of messages to
    stdout which are useful for debugging. The default value for verbose is 0.
    
      On success distances() returns a one-dimensional NumPy array which
    contains the Euclidean distances between all possible pairs of points in
    the input array. This array has form,
    
        dists = [1-to-2, 1-to-3, 1-to-4, ..., 2-to-3, 2-to-4, ...]
    
    where i-to-j is the Euclidean distance between points i and j. In general,
    for n input points, there will be (0.5 * n * (n - 1)) elements in this
    array. distances() carries out only the minimum number of calculations to
    find all pairwise distances over the set of input points.
    
      For the sake of both speed- and space-efficiency, the underlying memory
    for the array returned by distances() is a contiguous block of C doubles.
    On return, distances() transfers ownership of this memory to the returned
    NumPy arrray. As such this memory is freed by Python's garbage collector
    when all references to the returned NumPy array have vanished.
    
      The full signature for distances() is thus,
        
        distances(points[, threads = 1][, verbose = 0]) -> numpy.ndarray
    
    where points, threads and verbose may be specified positionally, by those
    keywords, or by a combination of the two.
    
      pairwise also provides the index() helper method. Given the size of a set
    of input points to distances(), n, and the indices of any two different
    points in that set, i and j, index() returns the index of the distance
    i-to-j in the corresponding array returned by distances(). The full
    signature for index() is thus,
    
        index(n, i, j) -> int
    
    where n, i and j are all integers, and n is greater than one. If i or j are
    negative, they index backwards from the end of the set of input points.
    That is, for i = 0 and j = -1, index() will return the index of the
    distance between the first and last points in the set of input points. i
    and j must not both refer to the same point. index() will raise a Python
    exception if any of these restrictions are violated.
    
      For further information on the inner-workings of pairwise please see
    source/pairwise.c, which is heavily commented.

