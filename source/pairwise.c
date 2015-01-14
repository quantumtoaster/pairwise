/*

    pairwise Version 1.0
    
    A Python extension module for the multi-threaded calculation of pairwise
    distances from a set of points in three-dimensional space.
    
    Copyright (C) 2015 Adam J. Marsden

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#define MIN(i, j) ((i < j) ? i : j)
#define MAX(i, j) ((i > j) ? i : j)

#define PW_VERSION "1.0"

#define PW_DISTANCES_DOCSTR "A function to calculate all the pairwise" \
                            " distances over a set of points in three-" \
                            "dimensional space. Only the minimum number of" \
                            " calculations is done, and these can be shared" \
                            " out amongst multiple threads. This function" \
                            " returns a new NumPy array of distances."

#define PW_INDEX_DOCSTR "A function to calculate the index into an array" \
                        " returned by pairwise.distances() at whose position" \
                        " the distance between any two different points with" \
                        " indices i and j will be found."

/*

    A structure for passing arguments to put_all_pairwise_distances() via
  pthread_create().
  
*/

typedef struct
put_all_pairwise_distances_arguments
{

    double* p;
    double* d;
    
    size_t n;
    size_t l;
    size_t u;
    
} papda_t;

/*

    A function to calculate all pairwise distances for a set of points in 
  three-dimensional space. The expected input is a contiguous array of doubles,
  p, with form,
  
    x1, y1, z1, x2, y2, z2, ...
  
  The output is a contiguous array of doubles, d, with form,
  
    d1to2, d1to3, ..., d2to3, ...
  
  The minimum number of calculations is done to find all pairwise distances.
  Given sane arguments the function will not fail, and will return (int)1 for
  success.
  
   Arguments are passed to this function via a papda_t structure. This allows
 the function to be called by pthread_create(), for running multiple pairwise
 distance calculations concurrently across multiple threads.
 
   Comment lines whose message begin "###" signal expansion points for the
 build.py script which accompanies this release.

*/

static int
put_all_pairwise_distances
(papda_t* args)
{

    /*
       Total number of input points.
    */
    register size_t n;
    
    /*
       Start, l (for lower limit), and end, u (for upper limit), points for a
       set of pairwise distance calculations.
    */
    register size_t l;
    register size_t u;
    
    /*
       Indices for the first, i, and second, j, points in the pair whose
       distance is to be calculated.
    */
    register size_t i;
    register size_t j;
    
    /*
       Indices for the position vectors of the first, a, and second, b, points
       in the pair whose distance is to be calculated.
    */
    register size_t a;
    register size_t b;
    
    /*
       Pointers to the input array of points, p, and the output array of
       distances, d.
    */
    register double* p;
    register double* d;
    
    /*
       x-, y-, and z-components of the position vector of the first point in
       the pair whose distance is to be calculated.
    */
    register double ax;
    register double ay;
    register double az;
    
    /*
       x-, y-, and z-components of the displacement vector from the first to
       the second point in the pair whose distance is to be calculated.
    */
    register double dx;
    register double dy;
    register double dz;
    
    //###UNROLLING_DECLARATIONS
    
    n = args->n;
    l = args->l;
    u = args->u;
    
    p = args->p;
    d = args->d;
    
    for (i = l; i < u; i ++) {
    
        a = i * 3;
            
        ax = *(p + a + 0);
        ay = *(p + a + 1);
        az = *(p + a + 2);
        
        j = i + 1;
        
        //###UNROLLING_ITERATIONS
    
        for (; j < n; j ++) {
        
            b = j * 3;
        
            dx = ax - *(p + b + 0);
            dy = ay - *(p + b + 1);
            dz = az - *(p + b + 2);
            
            *(d ++) = sqrt((dx * dx) + (dy * dy) + (dz * dz));
            
        }
        
    }
    
    return 1;

}

/*

    A function to populate papda_t structures with appropriately calculated
  variables such that a large set of pairwise calculations can be distributed
  roughly evenly across t separate smaller sets of pairwise calculations. It is
  expected that these calculations will then be carried out concurrently by t
  worker threads. Given sane arguments the function will not fail, and will
  return (int)1 for success.
  
    n is the total number of points which define the large set of pairwise
  calculations to be distributed across smaller sets. t is the desired number
  of these smaller sets of calculations. p is a pointer to the contiguous array
  which contains all n points, represented by triplets of doubles for each x-,
  y- and z-coordinate. As such the array pointed to by p should contain 3n
  elements. d is a pointer to a contiguous array to which all calculated
  pairwise distances, represented by doubles, will be written. As such the
  array pointed to by d should contain 0.5 * n * (n - 1) elements. a is a
  pointer to the contiguous array which contains all t papda_t structures, to
  which all calculated variables will be written.
  
    Note that this function does not check that the sizes of arrays pointed to
  by p, d or a are correct. Also note that this function does not modify the
  arrays pointed to by p or d; only the array pointed to by a is changed.

*/

static int
put_arguments_for_workers
(size_t n, size_t t, double* p, double* d, papda_t* a)
{
    
    /*
       Index for worker threads.
    */
    size_t i;
    
    /*
       Start, l (for lower limit), and end, u (for upper limit), points for a
       set of pairwise distance calculations.
    */
    size_t l;
    size_t u;
    
    /*
       Offset into the output distances array pointed to by d from which point
       the results of this set of pairwise distance calculations are to be
       written.
    */
    size_t d_offset;
    
    /*
       The ideal number of calculations to be done by each worker thread.
    */
    double c_per_t;
       
    double n_working;
    double x;
    
    /*
       Given [sum(n - k) from k = 1 to k = n] = 0.5 * n * (n - 1) = c, the
       total number of calculations to be done by all worker threads, then, 
       c_per_t = c / t, the ideal number of calculations to be done by each
       worker thread.
    */
    c_per_t = (0.5 * n * (n - 1)) / t;
    
    u = 0;
    
    n_working = n;
    
    for (i = 0; i < t; i ++) {
    
        if (!i) {
        
            /*
               The first thread output its calculated distances to the
               beginning of the output distances array, d.
            */
            d_offset = 0;
            
            /*
               The first threads starts its calculations at the first point.
            */
            l = 0;
        
        } else {
        
            /*
               Subsequent threads output their calculated distances to later
               parts of the output distances array, d.
               
               [sum(x - k) from k = 1 to k = y] = 0.5 * y * (2x - y - 1) with,
               
               x = (n - l), the number of first points whose pairwise distances
               haven't yet been calculated, and
                 
               y = (u - l), the number of first points whose pairwise distances
               this particular thread will calculate.
            */
            d_offset += 0.5 * (u - l) * ((2 * (n - l)) - (u - l) - 1);
            
            /*
               Subsequent threads start their calculations from where the
               previous thread left off.
            */
            l = u;
        
        }
        
        /* 
           An expression for the number of additional first points whose
           pairwise calculations this thread will do. This is not likely to be
           an integer.
        */
        x = 0.5 * (-sqrt((-8 * c_per_t) + (4 * n_working * n_working) - (4 * n_working) + 1) + (2 * n_working) - 1);
        
        /*
           The number of first points whose pairwise calculations are still to
           be done.
        */
        n_working -= x;
        
        /*
           Hence, the end point for this thread's set of calculations, rounded
           down to an integer. If this is the last thread then we guarantee it
           finishes off all remaining calculations.
        */
        u = (i == t - 1) ? n : u + x;
        
        /*
           Package the arguments calculated for this worker thread into the
           next papda_t.
        */
        (a + i)->p = p;
        (a + i)->d = d + d_offset;
        
        (a + i)->n = n;
        (a + i)->l = l;
        (a + i)->u = u;
        
    }
    
    return 1;
    
}

/*

    A function to find the number of available processors on any platform which
  provides /proc/cpuinfo. This value is written to the integer pointed to by t.
  If either /proc/cpuinfo cannot be read or contains no entries identified by
  the keyword "processor", this function returns (int)0 for failure. Otherwise
  this function returns (int)1 for success.

*/

static int
get_number_of_processors
(ssize_t* t)
{

    FILE* cpuinfo;
    
    char* line;
    
    size_t length;
    
    int count;
    
    cpuinfo = fopen("/proc/cpuinfo", "r");
    
    if (!cpuinfo) return 0;
    
    line = NULL;
    length = 0;
    count = 0;
    
    while (getline(&line, &length, cpuinfo) != -1) {
    
        if (strstr(line, "processor")) count ++;
    
    }
    
    fclose(cpuinfo);
    
    if (!count) return 0;
    
    *t = count;
    
    return 1;

}

/*

  pairwise.distances
    (
      points = [[x1, y1, z1], [x2, y2, z2], ...]
      [, threads = int]
      [, verbose = int]
    )

    A Python function to calculate all pairwise distances from points, a Python
  sequence of sequences of three numbers. On success this function returns a
  new NumPy array containing these distances, and whose underlying memory is an
  array of C doubles.
  
    Numbers contained in the input sequence of sequences represent the x, y and
  z coordinates respectively of the set of points across which all pairwise
  distances will be calculated. These are preliminarily copied into a new
  contiguous C array of doubles to speed up access during calculation. As such
  there is some small memory overhead to this function besides the large NumPy
  array of distances which it outputs. This additional overhead is negligible
  for all sizes of sets of input points for which the corresponding output
  array of distances can be stored in the memory available to modern day hosts.
  For example, for an input of 30k points, the memory overhead in copying all
  coordinates into a new intermediary array is 176 kB; the memory consumed by
  the output array of distances is 0.83 GB.
  
    During the writing of this function an attempt was made to speed it up by
  checking if the input Python object, containing the set of points across
  which all pairwise distances will be calculated, is already a NumPy array of
  the expected shape and with contiguous underlying memory from which point
  coordinates could be read directly. It was found that, in practice, these
  checks took slightly more time than the preliminary copy as described above
  for all reasonable sizes of sets of input points.
  
    By default all calculations are done by only one thread; a different number
  of concurrent calculation threads can be specified with threads. A special
  value of 0 for threads asks the function to use as many threads as are made
  available by the host, discovered by attempting to interrogate /proc/cpuinfo.
  The file system access inherent in threads = 0 has a significant time cost
  and should be avoided for a large number of repeat calls.
  
    By default nothing is printed to standard output; some messages which are
  useful for debugging are printed with verbose = 1.

*/

static PyObject*
py_distances
(PyObject* self, PyObject* values, PyObject* keys)
{

    /* 
       Python argument keywords.
    */
    static char* k[4] = {"points", "threads", "verbose", NULL};

    /*
       Pointers to the input, o, and output, r, Python objects; o should be a 
       sequence of sequences of three numbers, and r will be NumPy array.
    */
    PyObject* o; 
    PyObject* r;
    
    /*
       Number of threads to use for pairwise calculations, t, and level of
       verbosity, v.
    */
    ssize_t t;
    ssize_t v;
    
    /*
       Total number of input points.
    */
    Py_ssize_t n;
    
    /*
       Total number of pairwise calculations to be done.
    */
    size_t c;
    
    /*
       As c, but in an appropriate form to initialise a NumPy array.
    */
    npy_intp l[1];
    
    /*
       Pointers to the underlying data of the first, f, and second, s,
       dimensions of o.
    */
    register PyObject** f;
    register PyObject** s;
    
    /*
       Indices over f and s, respectively.
    */
    register size_t i;
    register size_t j;
    
    /*
       Pointers to contiguous arrays of input point coordinates, p, and output
       pairwise distances, d.
    */
    double* p;
    double* d;
    
    /*
       Pointer to an array of handles for each worker thread.
    */
    pthread_t* w;
    
    /*
       Pointer to an array of papda_t for each worker thread.
    */
    papda_t* a;
    
    /*
       Set default number of threads to use, t, to one.
    */
    t = 1;
    
    /*
       Set default verbosity level, v, to zero - report nothing.
    */
    v = 0;
    
    /*
       Get the input sequence of sequences of three floats, o, and, optionally,
       number of threads to use, t, and verbosity level, v. These argument
       values have keywords, "points", "threads" and "verbose", respectively.
       
       We accept long integers for t and v, and then ensure that these are
       positive later because the in-built exceptions raised by the PyArg
       function are not very helpful to the user; we can do better.
    */
    if (!PyArg_ParseTupleAndKeywords(values, keys, "O|ll:distances", k, &o, &t, &v)) return NULL;
    
    if (t < 0) {
        PyErr_SetString(PyExc_ValueError, "Supplied threads must be a positive integer.");
        return NULL;
    }
    
    if (v < 0) {
        PyErr_SetString(PyExc_ValueError, "Supplied verbose must be a positive integer.");
        return NULL;
    }
    
    /*
       If the user specifies the special value of zero for threads, attempt to
       read cpuinfo and use as many as processors are available. If then
       cpuinfo can't be read raise an exception.
    */
    if (!t) if (!get_number_of_processors(&t)) {
        PyErr_SetString(PyExc_RuntimeError, "Unable to detect the number of processors provided by the host.");
        return NULL;
    }
    
    if (!(o = PySequence_Fast(o, NULL))) {
        PyErr_SetString(PyExc_TypeError, "First dimension of supplied points does not support the sequence protocol.");
        return NULL;
    }
    
    /*
       Get the length of the input sequence of sequences of three floats.
    */
    if (!(n = PySequence_Length(o))) {
        PyErr_SetString(PyExc_ValueError, "First dimension of supplied points must not be an empty sequence.");
        return NULL;
    }
    
    /*
       Allocate memory for calculation input array, p. Free before return.
    */
    if (!(p = malloc(3 * n * sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for the calculation input array.");
        return NULL;
    }
    
    /*
       Expose, as f, the first dimension of the input sequence.
    */
    f = PySequence_Fast_ITEMS(o);
    
    for (i = 0; i < n; i ++) {
    
        if (!(f[i] = PySequence_Fast(f[i], NULL))) {
            PyErr_SetString(PyExc_TypeError, "Second dimension of supplied points does not support the sequence protocol.");
            return NULL;
        }
        
        /*
           Check that the second dimension of the input sequence has length
           three.
        */
        if (PySequence_Length(f[i]) != 3) {
            PyErr_SetString(PyExc_ValueError, "Second dimension of supplied points must have length three.");
            return NULL;
        }
        
        /*
           Expose, as s, the second dimension of the input sequence.
        */
        s = PySequence_Fast_ITEMS(f[i]);
        
        /*
           Check that the third dimension of the input the sequence is a
           number. If so, store that number as a double in the calculation
           input array, p.
        */
        for (j = 0; j < 3; j ++) {
            
            if (PyFloat_Check(s[j])) *(p ++) = PyFloat_AS_DOUBLE(s[j]);
            
            else if (PyInt_Check(s[j])) *(p ++) = (double)PyInt_AS_LONG(s[j]);
            
            else {
                PyErr_SetString(PyExc_ValueError, "Third dimension of supplied points must be a number.");
                return NULL;
            }
        
        }
    
    }
    
    /*
       Rewind pointer to calculation input array memory for later freeing.
    */
    p -= 3 * n;
    
    /*
       Get the total number of calculations to be done by the workers.
    */
    c = 0.5 * n * (n - 1);
    
    l[0] = c;
    
    /*
       Allocate memory for calculation output array, d. Will be returned as
       NumPy array.
    */
    if (!(d = malloc(c * sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for the calculation output array.");
        return NULL;
    }
    
    /*
       Allocate memory for worker arguments array, a. Free before return.
    */
    if (!(a = malloc(t * sizeof(papda_t)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for the worker arguments array.");
        return NULL;
    }
    
    /*
       Print some information if verbosity level is greater than zero.
    */
    if (v) printf("pairwise: points = %zu; calculations = %zu; threads = %zu\n", n, c, t);
    
    /*
       If we're only using one thread, we can forego both calculating bounds
       for each worker and using pthread.
    */
    if (t == 1) {
    
        a->p = p;
        a->d = d;
        
        a->n = n;
        a->l = 0;
        a->u = n;
        
        put_all_pairwise_distances(a);
    
    /*
       Otherwise we first calculate bounds for each worker such that the number
       of calculations each must do is roughly equal. Then we launch each
       worker using pthread, and block until they've all finished.
    */
    } else {

        /*
           Allocate memory for workers array, w. Free before return.
        */
        if (!(w = malloc(t * sizeof(pthread_t)))) {
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for the worker threads array.");
            return NULL;
        }
        
        /*
           Calculate and place arguments for workers in a.
        */
        put_arguments_for_workers(n, t, p, d, a);
        
        /*
           Launch workers.
        */
        for (i = 0; i < t; i ++) {
        
            pthread_create(w + i, NULL, (void*)put_all_pairwise_distances, a + i);
            
            if (v) printf("worker: i = %zu; n = %zu; l = %zu; u = %zu; p = %p; d = %p\n", 
                          i, (a + i)->n, (a + i)->l, (a + i)->u, (a + i)->p, (a + i)->d);
        
        }
        
        /*
           Block while all workers finish.
        */
        for (i = 0; i < t; i ++) pthread_join(w[i], NULL);
        
        /*
           Free memory for workers array, w.
        */
        free(w);
        
    }
    
    /*
       Free memory for worker arguments array, a.
    */
    free(a);
    
    /*
       Free memory for calculation input array, p.
    */
    free(p);
    
    /*
       Create a NumPy array, r, from the calculation output array, d. Nothing
       is copied.
    */
    r = PyArray_SimpleNewFromData(1, l, NPY_DOUBLE, d);
    
    /*
       Give ownership of the calculation output array, d, to r for later
       garbage collection in Python.
    */
    #if defined(NPY_ARRAY_OWNDATA)
    PyArray_ENABLEFLAGS((PyArrayObject*)r, NPY_ARRAY_OWNDATA);
    #else
    PyArray_UpdateFlags((PyArrayObject*)r, NPY_OWNDATA);
    #endif
    
    return r;

}

/*

  pairwise.index(n = int, i = int, j = int)

    A Python function which returns an integer suitable as an index into a
  distances array returned by pairwise.distances() (and incidentally also
  suitable for that returned by scipy.spatial.distance.pdist()). n is the
  number of points, for which there will be 0.5 * n * (n - 1) distances which
  are indexed by this function, and i and j are the indices of the points whose
  distance is indexed by the return value of this function.
  
    The value returned by this function is independent of the order of i and j.
  i and j may be negative in which case they index backwards from the end of
  the set of input points. For example, index(200, -1, -2) returns an index for
  the distance between points with indices 199 and 198. i and j must be non-
  equivalent indices. That is, i should not equal j, and, for example, given
  n = 5 and i = 4, j should not equal -1. Violations of these rules are caught
  and will raise a suitable Python exception.

*/

static PyObject*
py_index
(PyObject* self, PyObject* values, PyObject* keys)
{
    
    static char* k[4] = {"n", "i", "j", NULL};
    
    ssize_t n;
    ssize_t i;
    ssize_t j;
    
    ssize_t a;
    ssize_t b;
    ssize_t p;
    
    PyObject* r;
    
    if (!PyArg_ParseTupleAndKeywords(values, keys, "lll:index", k, &n, &i, &j)) return NULL;
    
    if (n < 2) {
        PyErr_SetString(PyExc_ValueError, "Supplied n must be a greater-than-one integer.");
        return NULL;
    }
    
    if (i < 0) i += n;
    
    if (i < 0) {
        PyErr_SetString(PyExc_IndexError, "Supplied i must be a valid index into an n-sized sequence.");
        return NULL;
    }
    
    if (j < 0) j += n;
    
    if (j < 0) {
        PyErr_SetString(PyExc_IndexError, "Supplied j must be a valid index into an n-sized sequence.");
        return NULL;
    }

    if (i == j) {
        PyErr_SetString(PyExc_IndexError, "Supplied i and j must not be equivalent indicies.");
        return NULL;
    }

    a = MIN(i, j);
    
    b = MAX(i, j);
    
    p = -a * (a - (2 * n) + 1);
    
    p >>= 1;
    
    p += b - a - 1;
    
    r = Py_BuildValue("n", p);
    
    return r;

}

static PyMethodDef
py_methods[] = {

	{"distances", (PyCFunction)py_distances, METH_VARARGS | METH_KEYWORDS, PW_DISTANCES_DOCSTR},
	{"index", (PyCFunction)py_index, METH_VARARGS | METH_KEYWORDS, PW_INDEX_DOCSTR},
	{NULL, NULL, 0, NULL}
	
};

PyMODINIT_FUNC
initpairwise(void)
{

    PyObject* m;
    
    m = Py_InitModule("pairwise", py_methods);
    
    PyModule_AddStringConstant(m, "__version__", PW_VERSION);
    
    import_array();

}

