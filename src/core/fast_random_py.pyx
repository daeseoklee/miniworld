from libc.stdlib cimport rand, RAND_MAX, srand
from libc.math cimport log, sqrt
import numpy as np
import cython
from time import time



cpdef void seed():
    srand(time())

cpdef void fix_seed():
    srand(100)

cpdef double uniform():
    cdef double r = rand()
    return r / RAND_MAX


cpdef int randint(int n,int m):
    cdef int k=m-n+1
    cdef double prob=1.0/k
    cdef int i=n
    cdef double u=uniform()
    while i<m:
        if u<prob:
            return i
        u-=prob
        i+=1
    return m

cpdef int randint_with_probs(int n,int m,double[:] probs):
    cdef int i=n
    cdef double u=uniform()
    while i<m:
        if u<probs[i-n]:
            return i
        u-=probs[i-n]
        i+=1
    return m
cpdef bint randbool():
    cdef double u=uniform()
    return u>0.5

@cython.boundscheck(False)
cpdef double[:] multiple_uniform(int n):
    cdef int i
    cdef double[:] result = np.empty(n, dtype='f8', order='C')
    for i in range(n):
        result[i] = uniform()
    return result

cpdef double gaussian():
    cdef double x1, x2, w

    w = 2.0
    while (w >= 1.0):
        x1 = 2.0 * uniform() - 1.0
        x2 = 2.0 * uniform() - 1.0
        w = x1 * x1 + x2 * x2

    w = ((-2.0 * log(w)) / w) ** 0.5
    return x1 * w

@cython.boundscheck(False)
cdef void assign_random_gaussian_pair(double[:] out, int assign_ix):
    cdef double x1, x2, w

    w = 2.0
    while (w >= 1.0):
        x1 = 2.0 * uniform() - 1.0
        x2 = 2.0 * uniform() - 1.0
        w = x1 * x1 + x2 * x2

    w = sqrt((-2.0 * log(w)) / w)
    out[assign_ix] = x1 * w
    out[assign_ix + 1] = x2 * 2



@cython.boundscheck(False)
cpdef double[:] multiple_gaussian(int n):
    cdef int i
    cdef double[:] result = np.empty(n, dtype='f8', order='C')
    for i in range(n // 2):  # Int division ensures trailing index if n is odd.
        assign_random_gaussian_pair(result, i * 2)
    if n % 2 == 1:
        result[n - 1] = gaussian()

    return result